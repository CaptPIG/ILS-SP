# ILS-SP for E-VRP-NL

A research-oriented C++20 implementation of iterated local search with set
partitioning (ILS-SP) for the electric vehicle routing problem with nonlinear
charging functions (E-VRP-NL).

The solver specializes the broader HEVRP-NL method of Wang et al. (2025) to the
Montoya et al. (2017) benchmark: one vehicle type, no customer demand, no
customer time windows, no charging-station waiting functions, and an objective
equal to driving time plus charging time.

> **Implementation status**
>
> This is an independent research implementation, not the authors' reference
> code. It passes the included contract and Gurobi integration tests, but it
> does not yet reproduce the solution quality and runtime reported in the TRC
> paper consistently, particularly on larger instances. Results from this
> repository should be identified as results of this implementation.

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

The implementation uses Gurobi instead of the CPLEX configuration reported by
Wang et al. (2025). Details that are not uniquely specified by the papers, such
as the charging path assigned to a gap created by a customer move, remain
explicit implementation choices. The default `paper` profile is the closest
implemented interpretation; all other search profiles are experimental.

## Repository layout

```text
include/ils_sp/    public headers
src/               library implementation and command-line executable
tests/             contract and Gurobi integration tests
data/evrp-2017/    minimal public benchmark subset used by the tests
```

Experiment runners, generated results, solver binaries, and article PDFs are
deliberately excluded from this source repository.

## Prerequisites

- CMake 3.24 or newer
- C++20-compatible compiler
- Gurobi C++ 10, 11, 12, or 13 with a valid license

CMake resolves Gurobi through `GUROBI_HOME`, falling back to
`/opt/gurobi-current` when available. Gurobi and its license are not included.

## Build

Development build with tests:

```bash
cmake --preset dev
cmake --build --preset dev --parallel 1
ctest --preset dev --parallel 1 --output-on-failure
```

Optimized build:

```bash
cmake --preset release
cmake --build --preset release --parallel 1
```

Set `-DILS_SP_BUILD_TESTS=OFF` at configure time when only the solver executable
is required.

## Quick start

```bash
build/release-make/ils-sp \
  --instance data/evrp-2017/tc0c10s2cf1.xml \
  --search-profile paper \
  --threads 1
```

Run `ils-sp --help` for the complete command-line interface. The paper-scale
defaults are 50,000 ILS iterations, a 5,000-route main pool, and an SP call
every 1,500 completed iterations. Per-SP engineering safeguards can be disabled
with `--sp-stall-limit 0 --sp-time-limit 0` when an uncapped run is required.

## Benchmark data

The repository includes only the four public Montoya instances required by the
test suite. The complete 120-instance E-VRP-NL benchmark is available from
[VRP-REP](https://www.vrp-rep.org/datasets/item/2016-0020.html) under dataset ID
`2016-0020`.

## License

The original source code and documentation in this repository are licensed
under the [BSD 3-Clause License](LICENSE). The benchmark XML files under
`data/evrp-2017/` are third-party research data and are not relicensed by this
project; see [data/README.md](data/README.md) for their source and attribution.

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
