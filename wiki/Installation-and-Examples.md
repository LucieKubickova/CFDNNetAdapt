# Installation and examples

This page describes how to install the dependencies used by CFDNNetAdapt and
how to run the examples included in this repository.

## What the repository contains

CFDNNetAdapt combines multi-objective evolutionary optimization with
deep-neural-network surrogate models. The repository contains:

- the framework implementation in [`src/`](../src/);
- local copies of the Platypus MOEA and pyrenn code in
  [`thirdParty/`](../thirdParty/);
- benchmark examples for ZDT and LZ test functions;
- a CFD-based converging-diverging shape-optimization example using OpenFOAM;
  and
- reference input, output, and indicator data in [`data/`](../data/).

The example scripts add `src/` and `thirdParty/` to `sys.path` themselves.
Run them from their example directory rather than from the repository root.

## Requirements

### Common requirements

- Linux or another Unix-like environment;
- Python 3.10 or newer is recommended (the repository contains Python 3.10
  bytecode and uses legacy, script-style Python modules);
- Python packages:
  - `numpy`
  - `scipy`
  - `dill`
  - `matplotlib` (required by post-processing and indicator scripts);
- Bash; and
- a C compiler and `make` if the OpenFOAM example is used.

Platypus and pyrenn are included locally under `thirdParty/`; they do not need
to be installed separately from PyPI.

### OpenFOAM requirements

The CFD example additionally requires:

- OpenFOAM 8, including `blockMesh`, `checkMesh`, `wmake`, and the
  `RunFunctions` helper;
- an OpenFOAM environment script available at
  `/opt/openfoam8/etc/bashrc`, as referenced by the supplied run scripts; and
- permission to compile the custom solver/library in
  `example/convDifShapeOptim/02_ofCodes/`.

If OpenFOAM is installed elsewhere, update the `source` line in
[`Allrun-test`](../example/convDifShapeOptim/10_baseCase/Allrun-test) and any
other run script that sources the OpenFOAM environment.

## Installation

Clone the repository and create an isolated Python environment:

```bash
git clone https://github.com/LucieKubickova/CFDNNetAdapt.git
cd CFDNNetAdapt
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy dill matplotlib
```

Confirm that the Python dependencies and local modules can be imported:

```bash
python - <<'PY'
import dill
import matplotlib
import numpy
import scipy
import sys

sys.path.insert(0, "src")
sys.path.insert(0, "thirdParty")
import CFDNNetAdaptV3
import pyrennModV3
print("CFDNNetAdapt dependencies imported successfully")
PY
```

No `requirements.txt` or packaging configuration is currently provided, so
the explicit `pip install` command above is the supported setup method.

## Running the benchmark examples

The benchmark cases evaluate analytical test functions, so they do not need
OpenFOAM. Each case contains:

- `Allrun.py`: full optimization run;
- `partialRun.py`: reduced/partial run;
- `testRun.py`: small smoke-test configuration;
- `computeIndicators.py`: indicator and comparison plots; and
- `postProcess.py`: post-processing of an optimization run.

Run a smoke test from the case directory. For example:

```bash
cd example/ZDTS/ZDT1_10
python testRun.py
```

Other ZDT cases are [`ZDT2_10`](../example/ZDTS/ZDT2_10),
[`ZDT3_10`](../example/ZDTS/ZDT3_10), [`ZDT4_02`](../example/ZDTS/ZDT4_02),
and [`ZDT6_02`](../example/ZDTS/ZDT6_02). LZ cases are under
[`example/LZS/`](../example/LZS/), with function families `LZ1` through `LZ5`,
`LZ7`, `LZ8`, and `LZ9`; each is provided with `_03` and `_10` variants.

The scripts use relative paths such as `00_prepData/` and write optimization
results to `01_algoRuns/`. Keep the supplied directory layout unchanged.
Runs can be computationally expensive; use `testRun.py` before starting
`Allrun.py`.

## Running the CFD shape-optimization example

The converging-diverging example is in
[`example/convDifShapeOptim/`](../example/convDifShapeOptim/). It couples the
optimizer to OpenFOAM case generation and simulation:

1. Enter the example directory and run the smoke-test configuration:

   ```bash
   cd example/convDifShapeOptim
   python testRun.py
   ```

2. For a full optimization, use `python Allrun.py`. Use
   `python partialRun.py` for the reduced configuration.

3. The generated OpenFOAM cases are placed in `ZZ_cases/`; logs and
   optimization output are written below the example directory.

The CFD scripts generate mesh and case files from
[`01_pyCodes/caseConstructor.py`](../example/convDifShapeOptim/01_pyCodes/caseConstructor.py).
The base OpenFOAM case is in [`10_baseCase/`](../example/convDifShapeOptim/10_baseCase/).
Compile the custom OpenFOAM code before running a case:

```bash
cd example/convDifShapeOptim/02_ofCodes
source /opt/openfoam8/etc/bashrc
./Allwmake
```

The CFD workflow is substantially more expensive than the analytical
benchmarks and requires a working OpenFOAM installation.

## Included data

The committed data directory contains 70 `.dat` files (about 4.8 MB):

### Benchmark results

Each of the five ZDT directories contains 13 files:

- `*_optSols.dat`: reference optimal solutions;
- `*_dnn_nondoms.dat`, `*_gpr_nondoms.dat`,
  `*_nsgaii_nondoms.dat`, and `*_socemo_nondoms.dat`: non-dominated
  solutions from the DNN, Gaussian-process, NSGA-II, and SOCEMO methods;
- matching `*_hvValues.dat` files: hypervolume histories; and
- matching `*_igdValues.dat` files: inverted-generational-distance histories.

The available ZDT datasets are [`ZDT1_10`](../data/ZDT1_10),
[`ZDT2_10`](../data/ZDT2_10), [`ZDT3_10`](../data/ZDT3_10),
[`ZDT4_02`](../data/ZDT4_02), and [`ZDT6_02`](../data/ZDT6_02).
The non-dominated solution files include CSV headers describing decision
variables (`x0`, `x1`, ...), objectives (`f1`, `f2`), and, where applicable,
surrogate predictions (`predf1`, `predf2`).

### CFD data

[`data/convDifShapeOptim/`](../data/convDifShapeOptim/) contains reference
results for the converging-diverging case:

- `default.dat`: default-case objective values;
- `baseline_common_all.dat` and `baseline_common_nondoms.dat`: baseline
  optimization results;
- `baseline_new_lastgens.dat`: final baseline generations; and
- `cfdnnetadapt_lastgens.dat`: final CFDNNetAdapt generations.

The original CFD sampling and scaled training data used by the example are
also kept under
[`example/convDifShapeOptim/00_prepCFDData/`](../example/convDifShapeOptim/00_prepCFDData/).
These files include all sampled solutions, non-dominated solutions, angular
scales, and the source CFD solutions.

## Troubleshooting

- **`ModuleNotFoundError`**: activate `.venv`, install the common packages,
  and run the script from its own example directory.
- **Missing `00_prepData` or `01_algoRuns` files**: restore the example's
  committed input directory and preserve the expected relative paths.
- **OpenFOAM command not found**: source the OpenFOAM 8 environment and check
  the path in the run script.
- **A CFD run stops during mesh checking**: inspect the generated case logs
  before increasing the optimization size; invalid geometry parameters are
  deliberately rejected by the example.
