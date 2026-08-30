#include "ils_sp/xml.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace ils_sp {
namespace {

struct XmlNode {
  std::string name;
  std::unordered_map<std::string, std::string> attributes;
  std::string text;
  std::vector<XmlNode> children;

  [[nodiscard]] const XmlNode* child(std::string_view child_name) const {
    const auto found = std::find_if(
        children.begin(), children.end(), [&](const XmlNode& candidate) {
          return candidate.name == child_name;
        });
    return found == children.end() ? nullptr : &*found;
  }

  [[nodiscard]] std::vector<const XmlNode*> children_named(
      std::string_view child_name) const {
    std::vector<const XmlNode*> result;
    for (const XmlNode& candidate : children) {
      if (candidate.name == child_name) {
        result.push_back(&candidate);
      }
    }
    return result;
  }
};

[[nodiscard]] std::string trim(std::string value) {
  const auto first = std::find_if_not(value.begin(), value.end(), [](char c) {
    return std::isspace(static_cast<unsigned char>(c)) != 0;
  });
  const auto last = std::find_if_not(value.rbegin(), value.rend(), [](char c) {
                      return std::isspace(static_cast<unsigned char>(c)) != 0;
                    }).base();
  if (first >= last) {
    return {};
  }
  return std::string(first, last);
}

[[nodiscard]] std::string decode_entities(std::string value) {
  const std::pair<std::string_view, std::string_view> replacements[] = {
      {"&amp;", "&"}, {"&lt;", "<"}, {"&gt;", ">"},
      {"&quot;", "\""}, {"&apos;", "'"}};
  for (const auto& [encoded, decoded] : replacements) {
    std::size_t position = 0;
    while ((position = value.find(encoded, position)) != std::string::npos) {
      value.replace(position, encoded.size(), decoded);
      position += decoded.size();
    }
  }
  return value;
}

class XmlParser {
 public:
  explicit XmlParser(std::string source) : source_(std::move(source)) {}

  [[nodiscard]] XmlNode parse_document() {
    skip_misc();
    XmlNode result = parse_element();
    skip_misc();
    if (position_ != source_.size()) {
      fail("trailing content after root element");
    }
    return result;
  }

 private:
  void fail(std::string_view message) const {
    throw std::invalid_argument("XML parse error at byte " +
                                std::to_string(position_) + ": " +
                                std::string(message));
  }

  [[nodiscard]] bool starts_with(std::string_view token) const {
    return std::string_view(source_).substr(position_).starts_with(token);
  }

  void skip_space() {
    while (position_ < source_.size() &&
           std::isspace(static_cast<unsigned char>(source_[position_])) != 0) {
      ++position_;
    }
  }

  void skip_until(std::string_view terminator) {
    const std::size_t end = source_.find(terminator, position_);
    if (end == std::string::npos) {
      fail("unterminated XML declaration or comment");
    }
    position_ = end + terminator.size();
  }

  void skip_misc() {
    while (true) {
      skip_space();
      if (starts_with("<?")) {
        skip_until("?>");
      } else if (starts_with("<!--")) {
        skip_until("-->");
      } else {
        return;
      }
    }
  }

  [[nodiscard]] std::string parse_name() {
    const std::size_t begin = position_;
    while (position_ < source_.size()) {
      const unsigned char value =
          static_cast<unsigned char>(source_[position_]);
      if (!(std::isalnum(value) != 0 || value == '_' || value == '-' ||
            value == ':' || value == '.')) {
        break;
      }
      ++position_;
    }
    if (begin == position_) {
      fail("expected XML name");
    }
    return source_.substr(begin, position_ - begin);
  }

  [[nodiscard]] std::string parse_quoted() {
    if (position_ >= source_.size() ||
        (source_[position_] != '\'' && source_[position_] != '\"')) {
      fail("expected quoted attribute value");
    }
    const char quote = source_[position_++];
    const std::size_t begin = position_;
    const std::size_t end = source_.find(quote, position_);
    if (end == std::string::npos) {
      fail("unterminated attribute value");
    }
    position_ = end + 1;
    return decode_entities(source_.substr(begin, end - begin));
  }

  [[nodiscard]] XmlNode parse_element() {
    if (position_ >= source_.size() || source_[position_] != '<' ||
        starts_with("</")) {
      fail("expected opening element");
    }
    ++position_;
    XmlNode node;
    node.name = parse_name();
    while (true) {
      skip_space();
      if (starts_with("/>")) {
        position_ += 2;
        return node;
      }
      if (starts_with(">")) {
        ++position_;
        break;
      }
      std::string key = parse_name();
      skip_space();
      if (position_ >= source_.size() || source_[position_] != '=') {
        fail("expected '=' after attribute name");
      }
      ++position_;
      skip_space();
      node.attributes.emplace(std::move(key), parse_quoted());
    }

    std::string accumulated_text;
    while (true) {
      if (position_ >= source_.size()) {
        fail("unterminated element");
      }
      if (starts_with("<!--")) {
        skip_until("-->");
        continue;
      }
      if (starts_with("</")) {
        position_ += 2;
        const std::string closing_name = parse_name();
        if (closing_name != node.name) {
          fail("mismatched closing element");
        }
        skip_space();
        if (position_ >= source_.size() || source_[position_] != '>') {
          fail("expected closing '>'");
        }
        ++position_;
        node.text = trim(decode_entities(std::move(accumulated_text)));
        return node;
      }
      if (source_[position_] == '<') {
        node.children.push_back(parse_element());
      } else {
        const std::size_t next = source_.find('<', position_);
        if (next == std::string::npos) {
          fail("unterminated text element");
        }
        accumulated_text.append(source_, position_, next - position_);
        position_ = next;
      }
    }
  }

  std::string source_;
  std::size_t position_{0};
};

[[nodiscard]] const XmlNode& required_child(const XmlNode& parent,
                                            std::string_view name,
                                            std::string_view context) {
  const XmlNode* result = parent.child(name);
  if (result == nullptr) {
    throw std::invalid_argument("missing " + std::string(context) + "/" +
                                std::string(name));
  }
  return *result;
}

[[nodiscard]] std::string required_text(const XmlNode& parent,
                                        std::string_view name,
                                        std::string_view context) {
  const std::string result = required_child(parent, name, context).text;
  if (result.empty()) {
    throw std::invalid_argument("empty " + std::string(context) + "/" +
                                std::string(name));
  }
  return result;
}

[[nodiscard]] double required_double(const XmlNode& parent,
                                     std::string_view name,
                                     std::string_view context) {
  const std::string value = required_text(parent, name, context);
  std::size_t consumed = 0;
  const double result = std::stod(value, &consumed);
  if (consumed != value.size() || !std::isfinite(result)) {
    throw std::invalid_argument("invalid number at " + std::string(context));
  }
  return result;
}

[[nodiscard]] int required_int_attribute(const XmlNode& node,
                                         std::string_view key,
                                         std::string_view context) {
  const auto found = node.attributes.find(std::string(key));
  if (found == node.attributes.end()) {
    throw std::invalid_argument("missing attribute " + std::string(key) +
                                " at " + std::string(context));
  }
  std::size_t consumed = 0;
  const int result = std::stoi(found->second, &consumed);
  if (consumed != found->second.size()) {
    throw std::invalid_argument("invalid integer attribute at " +
                                std::string(context));
  }
  return result;
}

[[nodiscard]] std::string required_attribute(const XmlNode& node,
                                             std::string_view key,
                                             std::string_view context) {
  const auto found = node.attributes.find(std::string(key));
  if (found == node.attributes.end() || found->second.empty()) {
    throw std::invalid_argument("missing attribute " + std::string(key) +
                                " at " + std::string(context));
  }
  return found->second;
}

}  // namespace

Instance load_montoya_instance(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::invalid_argument("cannot open instance " + path.string());
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  XmlNode root = XmlParser(buffer.str()).parse_document();
  if (root.name != "instance") {
    throw std::invalid_argument("unexpected XML root " + root.name);
  }

  const XmlNode& info = required_child(root, "info", "instance");
  const std::string dataset = required_text(info, "dataset", "info");
  const std::string name = required_text(info, "name", "info");
  if (name != path.stem().string()) {
    throw std::invalid_argument("instance name differs from filename");
  }

  const XmlNode& fleet = required_child(root, "fleet", "instance");
  const XmlNode& profile = required_child(fleet, "vehicle_profile", "fleet");
  const XmlNode& custom = required_child(profile, "custom", "vehicle_profile");
  const Vehicle vehicle{
      .speed_kmph = required_double(profile, "speed_factor", "vehicle_profile"),
      .consumption_wh_per_km =
          required_double(custom, "consumption_rate", "vehicle/custom"),
      .battery_capacity_wh =
          required_double(custom, "battery_capacity", "vehicle/custom"),
      .max_route_duration_h =
          required_double(profile, "max_travel_time", "vehicle_profile")};

  std::unordered_map<std::string, ChargingCurve> curves;
  const XmlNode& functions =
      required_child(custom, "charging_functions", "vehicle/custom");
  for (const XmlNode* function : functions.children_named("function")) {
    const std::string station_type =
        required_attribute(*function, "cs_type", "charging function");
    std::vector<ChargingPoint> points;
    for (const XmlNode* breakpoint : function->children_named("breakpoint")) {
      points.push_back(ChargingPoint{
          .time_h = required_double(*breakpoint, "charging_time", "breakpoint"),
          .energy_wh =
              required_double(*breakpoint, "battery_level", "breakpoint")});
    }
    if (!curves
             .emplace(station_type,
                      ChargingCurve(station_type, std::move(points)))
             .second) {
      throw std::invalid_argument("duplicate charging function " + station_type);
    }
  }

  std::unordered_map<int, double> service_by_customer;
  const XmlNode& requests = required_child(root, "requests", "instance");
  for (const XmlNode* request : requests.children_named("request")) {
    const int customer_id =
        required_int_attribute(*request, "node", "request");
    const double service =
        required_double(*request, "service_time", "request");
    if (!service_by_customer.emplace(customer_id, service).second) {
      throw std::invalid_argument("duplicate customer request");
    }
  }

  const XmlNode& network = required_child(root, "network", "instance");
  const XmlNode& xml_nodes = required_child(network, "nodes", "network");
  std::vector<Node> nodes;
  std::unordered_set<int> parsed_customers;
  for (const XmlNode* xml_node : xml_nodes.children_named("node")) {
    const int id = required_int_attribute(*xml_node, "id", "node");
    const int raw_type = required_int_attribute(*xml_node, "type", "node");
    if (raw_type < static_cast<int>(NodeKind::Depot) ||
        raw_type > static_cast<int>(NodeKind::Station)) {
      throw std::invalid_argument("unknown node type");
    }
    const NodeKind kind = static_cast<NodeKind>(raw_type);
    std::string station_type;
    double service_time = 0.0;
    if (kind == NodeKind::Station) {
      const XmlNode& node_custom =
          required_child(*xml_node, "custom", "station node");
      station_type = required_text(node_custom, "cs_type", "station/custom");
    } else if (kind == NodeKind::Customer) {
      const auto request = service_by_customer.find(id);
      if (request == service_by_customer.end()) {
        throw std::invalid_argument("customer has no request");
      }
      service_time = request->second;
      parsed_customers.insert(id);
    }
    nodes.push_back(Node{
        .id = id,
        .kind = kind,
        .x_km = required_double(*xml_node, "cx", "node"),
        .y_km = required_double(*xml_node, "cy", "node"),
        .service_time_h = service_time,
        .station_type = std::move(station_type)});
  }
  if (parsed_customers.size() != service_by_customer.size()) {
    throw std::invalid_argument("requests do not match customer nodes");
  }
  std::sort(nodes.begin(), nodes.end(),
            [](const Node& left, const Node& right) { return left.id < right.id; });
  return Instance(name, dataset, std::move(nodes), vehicle, std::move(curves));
}

}  // namespace ils_sp

