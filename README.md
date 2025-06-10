# Fluid-Slender Structure Interaction (FSI) Framework

## Synopsis

This repository contains a **Massively Parallel Fluid-Slender Structure Interaction Solver** designed for simulating interactions between many flexbile slender structures and unsteady fluid flow. The framework integrates:

- **FluTAS** (Fluid Transport Accelerated Solver)  
- An **Isogeometric geometrically exact Simo-Reissner (Cosserat) beam model**  
- An **immersed boundary method (IBM)** for coupling

It is optimized for **high-performance computing (HPC)** environments with support for **CPU architectures**. The tool is modular, extensible, and capable of handling large-scale coupled simulations in fluid dynamics and structural mechanics.

---

## Features

### a. FluTAS: Accelerated Fluid Solver
- Momentum equations advanced using 3rd-order Runge-Kutta (explicit diffusion)
- Pressure solved using an FFT-based direct solver
- Parallelized with:
  - **MPI** for CPUs

### b. Simo-Reissner Beam Model: Structural Solver
- Geometrically exact beam model (Cosserat rod theory)
- Spatial discretization via **arbitrary order NURBS-elements**
- Time integration using **2nd-order Newmark scheme**
- Supports static and dynamic simulations with gravity and external loads
- Pre- and post-processing routines for geometry and material properties.
- Parallelized with:
  - **MPI** for CPUs

### c. Immersed Boundary Method (IBM)
- A novel 2nd-order isogeometric/finite-difference IBM for coupling
- Supports arbitrary coarse structural meshes on fixed Eulerian grids
- Parallelized with:
  - **MPI** for CPUs

---

## Code Structure

The framework is organized into modular components:

- `apps/`– Application-specific main programs (`main__*.f90`)
- `src/` – Core fluid, structure, and coupling solvers
- `lib/` – Libraries for FFT and NURBS routines

---

## Requirements

- Fortran compiler with **MPI** and **OpenACC** support (e.g., GCC)
- **FFTW** library for Fourier transforms
- **NURBS library** (`igalib.f90`) for the structural solver

---

## Compilation and Usage

### Compilation
1. Install all prerequisites (see `REQ` file)
2. Use the provided `Makefile` to compile for your target architectur 

### Running Simulations
1. Prepare input files for fluid and structure solvers (`dns.in` and `param_fibm.f90`)
2. Run the executable with appropriate `mpirun` command

### Visualization
Use **ParaView** or similar tools to visualize outputs.

---

## Documentation

- `HOW_TO_USE` – Compilation and usage guide  
- `INFO_INPUT` – Description of input files  
- `INFO_VISU` – Visualization instructions  

---

## Outputs

### Fluid Solver:
- Velocity and pressure fields

### Structural Solver:
- Deformations, internal forces, stiffness matrices, inertial properties

### Coupled Interaction:
- Fluid-structure forces (two-way coupling)

---

## References

### FluTAS:
- P. Costa et al., *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows*, **Comput. Math. Appl.** (2018)  
- M. Crialesi-Esposito et al., *FluTAS: A GPU-accelerated finite difference code for multiphase flows* (2022)

### Cosserat Rod:
- V. Agrawal et al., *An efficient isogeometric/finite-difference immersed boundary method for the fluid–structure interactions of slender flexible structures*, **CMAME** (2024)

---

## Collaborate on This Project

We welcome contributions to enhance the framework!

- Have improvements? Submit a pull request

For questions or collaborations, contact **Dr. Vishal Agrawal**, **Prof. Artem Kulachenko**, and **Prof. Luca Brandt**.

---
