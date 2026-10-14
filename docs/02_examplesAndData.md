# Examples and data

Examples of CFDNNetAdapt usage contain:
- benchmark examples for ZDT and LZ test functions;
- a CFD-based shape-optimization example using OpenFOAM;
  and
- input, output, and indicator data discussed in the accompanying paper in [`data/`](../tree/main/data/).

The example scripts add `src/` and `thirdParty/` to `sys.path` themselves.
Run them from their example directory rather than from the repository root.

## Running the benchmark examples

The benchmark cases evaluate analytical test functions, so they do not need
OpenFOAM. Each case contains:

- `Allrun.py`: full optimization run (computationally expensive);
- `partialRun.py`: partial run with reduced number of iterations and seeds;
- `testRun.py`: small smoke-test configuration with minimal settings;
- `computeIndicators.py`: calculation of inverted generational distance and hypervolume indicators; and
- `postProcess.py`: visualization of optimization run results in the objective space.

All scripts write optimization results to `01_algoRuns/`, including log file, trained neural networks, and data from individual MOEA optimizations.

Example of running a smoke test for ZDT1 with 10 parameters from the case directory:

```bash
cd example/ZDTS/ZDT1_10
python3 testRun.py
```

Other ZDT cases are [`ZDT2_10`](../tree/main/example/ZDTS/ZDT2_10),
[`ZDT3_10`](../tree/main/example/ZDTS/ZDT3_10), [`ZDT4_02`](../tree/main/example/ZDTS/ZDT4_02),
and [`ZDT6_02`](../tree/main/example/ZDTS/ZDT6_02). LZ cases are under
[`example/LZS/`](../tree/main/example/LZS/), with functions `LZ1` through `LZ5`,
`LZ7`, `LZ8`, and `LZ9`; each is provided with 3 and 10 parameters.

## Running the CFD shape-optimization example

The CFD example is in
[`example/convDifShapeOptim/`](../tree/main/example/convDifShapeOptim/). It couples the
optimizer to OpenFOAM case generation and simulation.

The directory includes script for smote test `testRun.py`, reduced optimization run `partialRun.py` and full optimization run `Allrun.py`. In smoke test, the CFD simulations are only created and not evaluated. In the full optimization run, CFD simulations are evaluated in full and the whole optimization is computationally expensive.

All scripts write optimization results to `01_algoRuns/` and the generated OpenFOAM cases are placed in `ZZ_cases/`.

To create the OpenFOAM simulations, additional python scripts are placed in [`01_pyCodes`](../tree/main/example/convDifShapeOptim/01_pyCodes) and an OpenFOAM base case is saved in [`10_baseCase/`](../tree/main/example/convDifShapeOptim/10_baseCase/). To calculate cost function value, a custom OpenFOAM post-processing
utility has to be compiled:

```bash
cd example/convDifShapeOptim/02_ofCodes
source /opt/openfoam8/etc/bashrc
./Allwmake
```

## Included data

The committed data directory contains data from runs discussed
in the accompanying paper ([`DOI`](https://doi.org/10.1016/j.asoc.2026.114800)):

### Benchmark results

Each of the five ZDT directories contain:

- `*_optSols.dat`: reference optimal solutions;
- `*_dnn_nondoms.dat`, `*_gpr_nondoms.dat`,
  `*_nsgaii_nondoms.dat`, and `*_socemo_nondoms.dat`: non-dominated
  solutions from CFDNNetAdapt, CFDGPRAdapt, NSGA-II, and SOCEMO methods;
- matching `*_hvValues.dat` files: hypervolume histories; and
- matching `*_igdValues.dat` files: inverted-generational-distance histories.

The available ZDT datasets are [`ZDT1_10`](../tree/main/data/ZDT1_10),
[`ZDT2_10`](../tree/main/data/ZDT2_10), [`ZDT3_10`](../tree/main/data/ZDT3_10),
[`ZDT4_02`](../tree/main/data/ZDT4_02), and [`ZDT6_02`](../tree/main/data/ZDT6_02).

### CFD data

[`data/convDifShapeOptim/`](../tree/main/data/convDifShapeOptim/) contains
results for the CFD-based case:

- `default.dat`: default-geometry objective values;
- `baseline_common_all.dat` and `baseline_common_nondoms.dat`: baseline
  optimization results;
- `baseline_new_lastgens.dat`: results from next generations calculated with only MOEA; and
- `cfdnnetadapt_lastgens.dat`: final CFDNNetAdapt generations.

The original CFD sampling and scaled training data used by the example are
also kept under
[`example/convDifShapeOptim/00_prepCFDData/`](../tree/main/example/convDifShapeOptim/00_prepCFDData/).
These files include all sampled solutions, non-dominated solutions, allowed mesh deformation (min and max angles), and raw optimization results with un-evaluatable CFD solutions.