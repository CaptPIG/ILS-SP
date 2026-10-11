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
> This is an independent research prototype, not the authors' reference code.
> The current mainline is the frozen `439f_plus_exact_packed_DFS_mask` source:
> the 439f policy plus an exact packed-station DFS remaining-mask optimization.
> It does not consistently reproduce the solution quality and runtime reported
> in the TRC paper. Results must be identified as results of this implementation.
> The packed DFS change preserves ascending station enumeration and does not
> alter its floating-point arithmetic, pruning, dominance, descriptor order,
> or RNG usage. Focused equivalence checks are not a complete performance or
> paper-reproduction guarantee.

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
by the papers remain implementation choices. In particular, the retained typed
entry deletion and deficit-gap weighting policies are disclosed assumptions,
not established author behavior. The `paper` profile is the closest implemented
interpretation; it does not, by itself, select every option in the frozen
mainline run recipe below.

## Repository layout

```text
include/ils_sp/          public headers for the COPT mainline
src/                    byte-preserved core and command-line source
src/main_with_schedule.cpp  frontend with output-only retained-schedule export
SOURCE_PROVENANCE.json   frozen archive identity and per-source SHA-256 map
build_mainline.sh       archived GCC-compatible build recipe, adapted paths
```

Existing Gurobi source, test files, benchmark data, license, and citations remain
in the repository. The current build excludes the legacy Gurobi files and test
targets. This update does not claim that the historical test suite validates
the new mainline. Experiment runners, result logs, binaries, commercial SDKs,
licenses, and credentials are not included.

## Prerequisites

- Linux: the frozen source uses Linux filesystem APIs
- C++20 compiler; the frozen recipe used GCC 14.2.0
- COPT 8.0.7 C++ SDK and your own separately authorized license
- xxHash 0.8.3 runtime (`libxxhash.so.0`)
- CMake 3.24 or newer for the optional CMake build

Set `COPT_HOME` to the SDK directory containing `include/` and `lib/`.
The COPT SDK and license must be installed and configured separately; neither
is supplied by this repository. No setup command here downloads dependencies.

## Build

The direct build follows the frozen recipe (`-O3 -DNDEBUG -flto=2`) and creates
both frontends:

```bash
export COPT_HOME=/path/to/copt80
bash build_mainline.sh
build/mainline/ils-sp --help
build/mainline/ils-sp-with-schedule --help
```

`CXX` and `BUILD_DIR` can override the compiler and build directory.
All translation units use these four definitions:

```text
ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD=1
ILS_SP_FORWARD_ONLY_TWO_OPT=1
ILS_SP_TYPED_DELETE_ENTRY=1
ILS_SP_PACKED_DFS_MASK=1
```

Optional CMake builds use the same source files and definitions:

```bash
cmake --preset release -DCOPT_HOME="$COPT_HOME"
cmake --build --preset release --parallel 1
# Development build without LTO:
cmake --preset dev -DCOPT_HOME="$COPT_HOME"
cmake --build --preset dev --parallel 1
```

The CMake executables are under `build/release-make/` or `build/dev-make/`.
CMake does not run the legacy tests. A rebuilt executable has its own identity;
it must not be labeled with the historical binary SHA in the provenance file.

## Quick start

```bash
build/mainline/ils-sp \
  --instance data/evrp-2017/tc0c10s2cf1.xml \
  --search-profile paper --iterations 100 --seed 84001 --threads 1
```

This is a short CLI example. Use the complete recipe below when comparing the
frozen mainline. `ils-sp --help` lists the full interface; some available options
represent earlier experiments rather than the selected mainline configuration.

## Frozen mainline run recipe

Supply your benchmark XML file as `INSTANCE`. The documented N80 configuration
uses 50,000 iterations, Gamma 32, seed 84001, a 5,000-route main pool, SP and
backup intervals of 1,500, four SP-integer threads, and a 10-second solve-wide
active FRVCP wall budget:

```bash
INSTANCE=/path/to/benchmark.xml
export ILS_SP_CATALOG_RAM_BYTES=1073741824
export ILS_SP_CATALOG_BLOCK_BYTES=67108864
export ILS_SP_CATALOG_METADATA_BYTES=134217728
export ILS_SP_CATALOG_DISK_BYTES=4294967296

build/mainline/ils-sp-with-schedule \
  --instance "$INSTANCE" --iterations 50000 --seed 84001 --progress 1500 \
  --search-profile paper --promising-arcs 32 \
  --sp-stall-limit 0 --sp-time-limit 0 \
  --path-catalog complete --insertion-paths preserve-cs \
  --pool-identity coverage --continuation record --replace-path-attempts 1 \
  --partial-route-subset random --sp-frequency 1500 --backup-frequency 1500 \
  --pool-capacity 5000 --adaptive-sp-stall 0 \
  --route-opening-paths projected --initialization-paths physical \
  --preserve-swap-stations --sp-threads 4 --frvcp-wall-budget 10
```

The archived witness was additionally bounded externally to 600 CPU seconds,
720 wall seconds, and 4 GiB address space. Those process limits are not set by
the command above. The full N80 benchmark files are not added by this update.

To export the winning retained physical schedule, set
`ILS_SP_VALIDATION_JSON=/path/to/schedule.json` before using
`ils-sp-with-schedule`. This output-only frontend does not re-optimize the
solution. The output file is overwritten, so choose a new path for each run.

The optional disk tier is enabled only when `ILS_SP_CATALOG_DISK_DIR` is set to
an absolute, fresh directory path. The frozen implementation requires its
existing parent to be below `/workspace/`, rejects tmpfs/ramfs, and checks for
4,096-byte filesystem blocks and sufficient space. Unusable backing paths
fall back to the in-memory tier. Existing cache directories are never reused
or overwritten. Leave this variable unset for an in-memory run; merely setting
the disk-byte budget does not enable disk storage. This environment-specific
restriction is retained to preserve source exactness.

## Source identity

`SOURCE_PROVENANCE.json` maps every uploaded mainline source file to its path
and SHA-256 inside the frozen archive. Core and frontend bytes are unchanged;
only their repository paths and build metadata are adapted. The archive's
SHA-256 is:

```text
af49f1efffbebdf444e360253a2ee37911cbfe6cd34218a03ee56f2727243830
```

The mainline excludes the first-escape experiment and radix-sort changes.

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
