# ILS-SP for E-VRP-NL

A research-oriented C++20 implementation of iterated local search with set
partitioning (ILS-SP) for the electric vehicle routing problem with nonlinear
charging functions (E-VRP-NL).

The solver specializes the broader HEVRP-NL method of Wang et al. (2025) to the
Montoya et al. (2017) benchmark: one vehicle type, no customer demand, no
customer time windows, no charging-station waiting functions, and an objective
equal to driving time plus charging time.

**Implementation status (2026-10-10)**

This is an independent research implementation, not the authors' reference
code. The current core version completed one N40 instance/seed at the full
50,000-iteration budget; its N80 full-budget run is still in progress. These
checks do not establish full-benchmark quality, general speed improvements,
or reproduction of the paper's reported runtime. Results should be identified
as results of this implementation.

## Algorithm provenance

The target method is the ILS-SP algorithm presented in the 2025 TRC paper. Its
search design builds on the Hybrid LNS framework from the 2023 EJOR paper. The
table below records the primary source used for each implemented component.

| Implemented component | Primary source |
| --- | --- |
| ILS-SP control flow and Route Generator/Route Assembler alternation | Wang et al. (2025), Algorithm 1 |
| Generalized objective, feasible/infeasible exploration, and adaptive penalty | Wang et al. (2025), Section 5.1 and Algorithm 2; Wang and Zhao (2023), Sections 4.1-4.2 |
| Fast route evaluation with nonlinear charging | Wang et al. (2025), Section 5.1.2 and Algorithm 3 |
| Initial solution and destroy-and-repair perturbation | Wang and Zhao (2023), Sections 4.3-4.4; adapted in Wang et al. (2025), Sections 5.1.3-5.1.4 |
| VND with Swap, Relocate, 2-OPT\*, and ReplacePath | Wang and Zhao (2023), Section 4.5; adapted in Wang et al. (2025), Section 5.1.5 |
| Promising-arc granular neighborhoods and customer correlation | Vidal et al. (2013); adopted in Wang and Zhao (2023), Section 4.6, and Wang et al. (2025), Sections 5.1.4-5.1.5 |
| Route pool, set partitioning, LP reduced costs, and FRVCP | Wang et al. (2025), Section 5.2 and Algorithm 5; route-pool/SP predecessor in Wang and Zhao (2023), Algorithm 1 |
| E-VRP-NL benchmark and piecewise-linear charging data | Montoya et al. (2017) |

The implementation uses COPT instead of the CPLEX configuration reported by
Wang et al. (2025). Details that are not uniquely specified by the papers, such
as the charging path assigned to a gap created by a customer move, remain
explicit implementation choices. The default `paper` profile is the closest
implemented interpretation; all other search profiles are experimental.

## Repository layout

```text
include/ils_sp/    public headers
src/              core library and command-line executable
CMakeLists.txt    COPT/xxHash build configuration
CMakePresets.json development and release presets
```

This snapshot contains only the current core source. Tests, benchmark data,
experiment runners, generated results, binaries, caches, and solver SDKs are
excluded.

## Prerequisites

- Linux, with a C++20-compatible compiler and Make
- CMake 3.24 or newer
- COPT 8.0.7 C++ SDK and a valid COPT license, installed separately
- A separately installed xxHash shared library; version 0.8.3 is required for
  the optional indexed disk-catalogue backend

Install and configure COPT using the [official installation guide](https://guide.coap.online/copt/en-doc/install.html).
Set `COPT_HOME` to its installation directory. CMake searches for `coptcpp.h`,
`libcopt_cpp`, and `libcopt` under that directory. The official SDK's
`include/coptcpp_inc` directory is also required. Configure the license outside
this repository and ensure the COPT shared libraries are available at runtime.
COPT's SDK and license are not covered by this project's BSD license.

The xxHash API declarations in this repository link to the installed library;
no hash implementation is bundled. See [xxHash 0.8.3](https://github.com/Cyan4973/xxHash/tree/v0.8.3).
If CMake cannot find it, pass `-DILS_SP_XXHASH_LIBRARY=...` with the full path
to the installed shared library.

## Build

From the repository root, after setting `COPT_HOME`:

```bash
cmake --preset release -DCOPT_HOME="$COPT_HOME"
cmake --build --preset release --parallel 1
```

Use the `dev` preset for a debug build. Both presets build the core library and
`ils-sp` executable only. The source retains its existing Linux-specific
filesystem implementation; portability to Windows or macOS is not claimed.

## Run

Supply a Montoya-format XML instance obtained separately. For example, save
one as `input.xml` and run the explicit current configuration:

```bash
build/release-make/ils-sp \
  --instance input.xml \
  --search-profile paper \
  --seed 84001 \
  --iterations 50000 \
  --threads 1 \
  --sp-threads 4 \
  --promising-arcs 40 \
  --sp-stall-limit 0 \
  --sp-time-limit 0 \
  --progress 1500 \
  --path-catalog complete \
  --insertion-paths completed \
  --pool-identity coverage \
  --continuation record \
  --replace-path-attempts 1 \
  --partial-route-subset random \
  --sp-frequency 1500 \
  --backup-frequency 1500 \
  --pool-capacity 5000 \
  --adaptive-sp-stall 0 \
  --route-opening-paths projected \
  --initialization-paths physical \
  --preserve-swap-stations
```

`--threads 1 --sp-threads 4` uses four threads only for the integer SP model;
LP and FRVCP retain one thread. Omitting options can select different CLI
defaults. Run `build/release-make/ils-sp --help` for the full interface.

The indexed disk-catalogue backend is opt-in via `ILS_SP_CATALOG_DISK_DIR`.
Leave that variable unset for ordinary use. The unchanged backend requires
xxHash 0.8.3 and a new directory whose existing parent resolves beneath
`/workspace/`, on a non-memory filesystem with 4096-byte blocks and sufficient
free space. These guards are retained from the frozen core version.

## Benchmark data

Benchmark XML files are not included. The Montoya E-VRP-NL benchmark is
available from [VRP-REP](https://www.vrp-rep.org/datasets/item/2016-0020.html),
dataset ID `2016-0020`. Refer to that source for attribution and terms.

## License

The original source code and documentation in this repository are licensed
under the [BSD 3-Clause License](LICENSE). COPT and xxHash are separate
dependencies governed by their respective licenses.

## References

1. W. Wang, Y. Adulyasak, J.-F. Cordeau, and G. He, “The
   Heterogeneous-Fleet Electric Vehicle Routing Problem with Nonlinear Charging
   Functions,” *Transportation Research Part C: Emerging Technologies*, 170,
   104932, 2025. <https://doi.org/10.1016/j.trc.2024.104932>
2. W. Wang and J. Zhao, “Partial linear recharging strategy for the electric
   fleet size and mix vehicle routing problem with time windows and recharging
   stations,” *European Journal of Operational Research*, 308(2), 929-948,
   2023. <https://doi.org/10.1016/j.ejor.2022.12.011>
3. T. Vidal, T. G. Crainic, M. Gendreau, and C. Prins, “A hybrid genetic
   algorithm with adaptive diversity management for a large class of vehicle
   routing problems with time-windows,” *Computers & Operations Research*,
   40(1), 475-489, 2013. <https://doi.org/10.1016/j.cor.2012.07.018>
4. A. Montoya, C. Guéret, J. E. Mendoza, and J. G. Villegas, “The electric
   vehicle routing problem with nonlinear charging function,” *Transportation
   Research Part B: Methodological*, 103, 87-110, 2017.
   <https://doi.org/10.1016/j.trb.2017.02.004>
