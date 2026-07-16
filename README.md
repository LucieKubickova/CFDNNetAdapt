# CFDNNetAdapt
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.729050990.svg)](https://doi.org/10.5281/zenodo.21158857)

CFDNNetAdapt is an adaptive CFD–DNN optimization framework for CFD-based shape optimization. The algorithm combines computational fluid dynamics (CFD) simulations with multi-objective evolutionary optimization (MOEA) and directly integrated deep neural network (DNN) surrogate models. The DNN architecture is identified automatically during the optimization process and is used to accelerate the computationally expensive mid-to-late optimization stages.

This repository accompanies the paper:

> **Lucie Kubíčková, Ondřej Gebouský, Jan Haidl, Martin Isoz**  
> *Accelerating shape optimization by deep neural networks with on-the-fly determined architecture*  
> Applied Soft Computing, Volume 192, 2026, 114800  
> DOI: https://doi.org/10.1016/j.asoc.2026.114800

## Repository Description

The repository contains the implementation of the CFDNNetAdapt methodology used for accelerating multi-objective CFD-based shape optimization. The framework combines:

- Multi-objective evolutionary algorithms (MOEA)
- Adaptive deep neural network surrogate models
- Automatic DNN architecture selection
- CFD-based objective-function evaluation
- Benchmark optimization cases
- Single-phase ejector shape optimization workflows

## Third-Party Software

### MOEA

D. Hadka, *Platypus: A Free and Open Source Python Library for Multiobjective Optimization*, 2020.

https://github.com/Project-Platypus/Platypus

### Neural Networks

D. Atabay, *pyrenn: A recurrent neural network toolbox for Python and MATLAB*, Technische Universität München, 2018.

https://pyrenn.readthedocs.io/en/latest/

## Citation

If you use this software, please cite:

```text
Kubíčková, L., Gebouský, O., Haidl, J., & Isoz, M. (2026).
Accelerating shape optimization by deep neural networks with on-the-fly determined architecture.
Applied Soft Computing, 192, 114800.
https://doi.org/10.1016/j.asoc.2026.114800
