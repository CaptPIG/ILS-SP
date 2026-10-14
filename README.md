# ILS-SP for E-VRP-NL

A research-oriented C++20 implementation of iterated local search with set
partitioning (ILS-SP) for the electric vehicle routing problem with nonlinear
charging functions (E-VRP-NL).

## Frozen R4 mainline

The current mainline publishes the frozen R4 source variant
`n40_dual_start_and_rc_ub0_v2`. The 27 core and command-line source files are
byte-preserved from the frozen input; only repository paths and build metadata
are adapted. This is an independent research prototype, not the authors'
reference implementation or a claim of reproducing the paper's results.

R4 integrates reduced-cost fixing for the N > 10 paper-profile unlimited-SP
entry using the existing full LP. The full route pool and original LP feedback
remain intact; fixing changes only the current MIP upper bounds. N <= 10 skips
this prepass. The frozen configuration disables typed entry deletion and retains
physical reverse-empty-child handling, forward-only two-opt, and packed DFS masks.
No later two-opt profile-reuse experiment is included.

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

The current mainline uses COPT 8.0.7 instead of the CPLEX configuration
reported by Wang et al. (2025). Earlier Gurobi files are retained for historical
context but are not part of the current build. Details not uniquely specified
by the papers remain implementation choices. In particular, the retained deficit-gap weighting policy is a disclosed assumption,
not established author behavior. R4 disables typed entry deletion. The `paper` profile is the closest implemented
interpretation; it does not, by itself, select every option in the frozen
mainline run recipe below.


## Repository layout

- `include/ils_sp/`: frozen public headers
- `src/`: frozen core and command-line source, alongside retained legacy files
- `src/main_with_schedule.cpp`: frontend with retained-schedule output
- `SOURCE_PROVENANCE.json`: frozen identity and per-source SHA-256 map
- `build_mainline.sh`: frozen GCC-compatible build flags with repository paths

Existing Gurobi source, test files, benchmark data, license, and citations are
retained. The current build excludes legacy Gurobi files and test targets.
This publication adds no test artifacts, experiment runners, result logs,
binaries, commercial SDKs, license keys, or local recovery manifests.

## Prerequisites

- Linux and a C++20 compiler; the frozen build used GNU C++ 14.2.0
- COPT 8.0.7 C++ SDK and your own separately authorized license
- xxHash runtime (`libxxhash.so.0`)
- CMake 3.24 or newer for the optional CMake build

Install dependencies separately and set `COPT_HOME` to the SDK directory
containing `include/` and `lib/`. No command below downloads dependencies.

## Build

The direct recipe uses `-O3 -DNDEBUG -flto=1` and builds both frontends:

```bash
export COPT_HOME=/path/to/copt80
bash build_mainline.sh
build/mainline/ils-sp --help
build/mainline/ils-sp-with-schedule --help
```

`CXX` and `BUILD_DIR` may override the compiler and output directory.
All translation units use the frozen R4 definitions:

```text
ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD=1
ILS_SP_FORWARD_ONLY_TWO_OPT=1
ILS_SP_TYPED_DELETE_ENTRY=0
ILS_SP_PACKED_DFS_MASK=1
```

The existing CMake release preset uses the same source and definitions:

```bash
cmake --preset release -DCOPT_HOME="$COPT_HOME"
cmake --build --preset release --parallel 1
```

Executables are written under `build/release-make/`. The `dev` preset disables
LTO and enables a Debug build; it is not the frozen optimized recipe.
Build metadata was adapted for this publication; no new build or tests were
run as part of the upload. A rebuilt executable has its own identity.

## Usage

```bash
build/mainline/ils-sp \
  --instance data/evrp-2017/tc0c10s2cf1.xml \
  --search-profile paper --iterations 100 --seed 84001 --threads 1
```

This is a short CLI example, not a benchmark recipe. Use `--help` to inspect
the complete interface and set all relevant options for your experiment.

To export the retained physical schedule, set
`ILS_SP_VALIDATION_JSON=/path/to/schedule.json` and run
`ils-sp-with-schedule`. This output frontend does not re-optimize the solution.
The output file is overwritten, so choose a fresh path for each run.

The optional disk tier is enabled only when `ILS_SP_CATALOG_DISK_DIR` is set
to an absolute fresh directory path. The frozen implementation requires its
existing parent to be below `/workspace/`, rejects tmpfs/ramfs, and checks for
4,096-byte filesystem blocks and sufficient space. Unusable backing paths
fall back to the in-memory tier. Leave the variable unset for an in-memory run.
This environment-specific restriction is preserved in the original source.

## Source identity

The frozen 44-file input identity has SHA-256:

```text
a43f547aebfc8841eef83b69fb3d5e4a7842867b9347f4ea5f5d689a23272633
```

`SOURCE_PROVENANCE.json` identifies the 27 published core/frontend source files,
their relative frozen-input paths, and byte-level hashes. The remaining frozen
input files are recovery metadata, tests, and patches that are not published.
The source identity is not an archive hash or an executable hash.

## Benchmark data

The repository retains the four public Montoya instances from the original
test subset. The complete 120-instance E-VRP-NL benchmark is available from
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
