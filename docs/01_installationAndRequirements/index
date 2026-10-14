# Installation and requirements

This page describes how to install the dependencies used by CFDNNetAdapt.

## What the repository contains

CFDNNetAdapt combines multi-objective evolutionary optimization with
deep-neural-network surrogate models. The repository contains:

- the framework implementation in [`src/`](../tree/main/src/);
- local copies of the Platypus MOEA and pyrenn code in
  [`thirdParty/`](../tree/main/thirdParty/).

## Requirements

### Common requirements

- Python 3.10 or newer;
- Python packages:
  - `numpy`
  - `scipy`
  - `dill`
  - `matplotlib` (optional visualization of results);

Platypus and pyrenn are included locally under `thirdParty/`; they do not need
to be installed separately from PyPI.

### OpenFOAM requirements

The CFD example additionally requires:

- OpenFOAM v8;
- an OpenFOAM environment script available at
  `/opt/openfoam8/etc/bashrc`, as referenced by the supplied run scripts.

If OpenFOAM is installed elsewhere, update the `source` line in
[`Allrun-serialAll`](../tree/main/example/convDifShapeOptim/10_baseCase/Allrun-serialAll) and [`Allrun-test`](../tree/main/example/convDifShapeOptim/10_baseCase/Allrun-test).

## Installation

Clone the repository and create an isolated Python environment:

```bash
git clone https://github.com/LucieKubickova/CFDNNetAdapt.git
cd CFDNNetAdapt
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy scipy dill matplotlib
```

Confirm that the Python dependencies and local modules can be imported:

```bash
python3 - <<'PY'
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