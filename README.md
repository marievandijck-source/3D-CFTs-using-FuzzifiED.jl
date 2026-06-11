# 3D CFTs using FuzzifiED.jl

This repository contains the code used for the second part of my master's thesis at Ghent University, on the topic of three-dimensional conformal field theories (CFTs).

## Overview

The model under investigation is the **O(3) Wilson-Fisher CFT**, studied using the quantum rotor model on the fuzzy sphere. The numerical calculations were performed using [FuzzifiED.jl](https://github.com/FuzzifiED/FuzzifiED.jl).

The approach follows the methods described in:
> Zheng Zhou et al., *FuzzifiED: A Julia Package for Fuzzy Sphere Calculations*, [arXiv:2510.09755](https://arxiv.org/abs/2510.09755)

Two main quantities were calculated:

- **Scaling dimensions** of CFT operators, obtained via both Exact Diagonalization (ED) and the Density Matrix Renormalization Group (DMRG)
- **OPE coefficients**, extracted from the eigenstates

More details on each calculation can be found in the associated subfolders.

## Repository Structure

| Folder | Contents |
|---|---|
| `Scaling Dimensions` | ED and DMRG code for scaling dimensions, figure scripts, and output figures |
| `OPE` | Code for extracting OPE coefficients from the eigenstates |

## Dependencies

- `JLD2.jl`, `Plots.jl`, `LaTeXStrings.jl`, `LinearAlgebra.jl`, `ITensors.jl`, `ITensorMPS.jl`
- [FuzzifiED.jl](https://github.com/FuzzifiED/FuzzifiED.jl)
