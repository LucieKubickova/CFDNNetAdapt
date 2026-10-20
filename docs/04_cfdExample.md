# CFD shape-optimization example

This section describes the CFD-based shape-optimization example from
[`example/convDifShapeOptim/`](../tree/main/example/convDifShapeOptim/). It couples
CFDNNetAdapt to OpenFOAM simulations of a single-phase ejector.

## Physical problem and objectives

The case models a single-phase ejector with an inlet, suction, converging part,
mixing tube, and diffuser. The geometry is controlled by 14 parameters
consisting of:

- $$x$$ and $$y$$ coordinates of two control points defining the shape of the
  converging part;
- $$x$$ and $$y$$ coordinates of four control points defining the shape of the
  diffuser;
- length of the converging part ($$\ell_\mathrm{C}$$); and
- length of the diffuser ($$\ell_\mathrm{D}$$).

Candidate geometries are rejected when control-point x-coordinates are not
strictly increasing or when

$$
\ell_\mathrm{C}+\ell_\mathrm{D}\leq\ell_\mathrm{max},
$$

where $$\ell_\mathrm{max}$$ is a predefined maximal value.

The optimized objectives include the ejector total length

$$
\ell_\mathrm{t}=\ell_\mathrm{C}+\ell_\mathrm{D},
$$

and the ejector energy efficiency calculated from flow rates $$\dot{Q}$$ and
pressure values $$p$$ as

$$
e_\mathrm{eff} =
\frac{Q_\mathrm{suction}(p_\mathrm{D}-p_\mathrm{suction})}
{Q_\mathrm{inlet}(p_\mathrm{inlet}-p_\mathrm{D})}.
$$

The optimization maximizes $$e_\mathrm{eff}$$ and minimizes $$\ell_\mathrm{t}$$.

The energy efficiency is calculated by the custom OpenFOAM utility
[`calcEnergyEfficiencyV1.C`](../tree/main/example/convDifShapeOptim/02_ofCodes/calcEnergyEfficiencyV1.C).
The utility has to be compiled before running the optimization:

```bash
cd example/convDifShapeOptim/02_ofCodes
source /opt/openfoam8/etc/bashrc
./Allwmake
```

## Running the optimization

Full optimization run:

```bash
cd example/convDifShapeOptim
python Allrun.py
```

## `Allrun.py` settings

[`Allrun.py`](../tree/main/example/convDifShapeOptim/Allrun.py) configures a
CFDNNetAdaptV3 run as follows:

- `nPars = 14` -- number of decision variables
- `nObjs = 2` -- number of objectives
- `nOuts = 1` -- number of surrogate outputs
- `nSam = 2000` -- initial sample count
- `deltaNSams = [500]` -- sample count increments over iterations
- `nNN = 10` -- number of neural networks generated in each iteration
- `minN = 2` -- minimum number of neurons in a hidden layer
- `maxN = 20` -- maximum number of neurons in a hidden layer
- `nHidLay = 3` -- number of hidden layers
- `dRN = 0` -- decrement of variability in hidden layer size
- `tol = 5e-2` -- optimization tolerance
- `iMax = 200` -- maximum number of iterations
- `nComps = 8` -- number points chosen for verification
- `nSeeds = 5` -- number of seeds for DNN training
- `trainPro = 75` -- training percentage split
- `valPro = 15` -- validation percentage split
- `testPro = 10` -- testing percentage split
- `kMax = 10000` -- maximum number of iterations in DNN training
- `rEStop = 1e-5` -- relative error for stopping DNN training
- `popSize = 500` -- population size of MOEA algorithm
- `offSize = 500` -- offspring size of MOEA algorithm
- `nGens = 30` -- number of MOEA generations

## Geometry generation and OpenFOAM execution

Scripts for generation of CFD simulations are stored in
[`01_pyCodes/`](../tree/main/example/convDifShapeOptim/01_pyCodes/).
[`caseConstructor.py`](../tree/main/example/convDifShapeOptim/01_pyCodes/caseConstructor.py)
copies the base case, generates `blockMeshDict`, and copies the current
geometry scripts into the generated case. The base case
[`10_baseCase/`](../tree/main/example/convDifShapeOptim/10_baseCase/) supplies the
boundary conditions, turbulence model, solver settings, and field files.

The CFD simulation is run by executing `Allrun-serialAll`, which:

- sources OpenFOAM v8;
- generates the mesh by running `blockMesh`;
- checks the mesh by running `checkMesh`;
- evaluates the simulation via `simpleFoam`; and
- calculates the energy efficiency with `calcEnergyEfficiencyV1`.

Simulations with failed mesh generation, missing computation, failed residual
checks, or missing/invalid energy-efficiency output are rejected and not
included in the optimization.

## Outputs

`Allrun.py` saves data to `01_algoRuns` in a dedicated `run_*` folder,
including a `log.out` file and folders for each optimization step `step_*`.

The `log.out` file stores data about:

- used settings;
- architectures of all tested DNNs;
- architecture of the best-found DNN;
- points chosen for verification with DNN-predicted values and values
  calculated by CFD; and
- verification errors.

Each `step_*` folder includes:

- weights and metadata of each trained DNN; and
- MOEA results from each DNN-based optimization.

Individual CFD simulations are not stored by default, but the simulation
folder is cleared after each iteration. To save the CFD simulations, set
`cleanSims = False` in `Allrun.py`. However, when running more iterations of
CFDNNetAdapt, the saved simulations may be overwritten.

## Visualization of results

[`postProcess.py`](../tree/main/example/convDifShapeOptim/postProcess.py) plots results
in the objective space, including:

- the NSGA-II sample data; and
- the solutions predicted by the best DNNs found in each iteration.

The resulting figure `objSpacePlot.png` is saved to the run folder.

![Figure generated by `postProcess.py`](https://github.com/LucieKubickova/CFDNNetAdapt/tree/main/Images/objSpacePlotCFD.png)