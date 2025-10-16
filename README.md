# CurvatureEngine

## Overview
This project is a high-performance computational framework for solving geodesic equations and curvature tensors in General Relativity. It includes tools for computing Christoffel symbols, Ricci tensors, Riemann curvature tensors, and the evolution of spacetime using the ADM/BSSN formalism. The code supports multiple metrics, including Schwarzschild and Minkowski, and provides efficient numerical solvers. In the future, this solver will Use the entire Tensorium_lib and the initial condition and generated code and optimizations by our MLIR compiler (Tensorium_MLIR).
## Features
- **Geodesic Integrator:** Computes null and timelike geodesics in curved spacetimes.
- **Curvature Tensor Computation in analytical formalism:** Computes Christoffel symbols, Riemann, Ricci, and Einstein tensors.
- **BSSN Formalism Solver:** Implements evolution equations for numerical relativity using a 3+1 decomposition.
- **Constraint Equations:** Enforces Hamiltonian and momentum constraints to ensure numerical stability.
- **Support for Multiple Metrics:** Includes Schwarzschild, Minkowski, and Kerr metrics, with extensibility to other spacetimes.
- **Parallelized Computation:** Optimized for multi-threading and SIMD vectorization to accelerate computations (Parallel computing is not ready atm so monothread atm...).
- **High-Precision Numerical Methods:** Uses high-order finite difference schemes for derivatives and tensor calculations.

## Compilation
This project is optimized for **Clang++** with advanced compiler flags for performance:
```sh
make
```

## Usage
Run the simulation with:
```sh
./CurvatureEngine <options>
```
Data will be stored in CSV and VTK (for extrinsic curvature tensor components and Geodesics)format for further analysis.
The plots are managed in the scripts (all python and for k_to_vtk you need paraview and pvpython)
## Future Work (TODO)
- Aarch64 support on Macos 
- Unit tests and debug mode.
- Class template to integrate your own problems
- Optimize for parallel computation and SIMDs.
- Implement relativistic fluid dynamics template.
- Implement GPU acceleration for faster computations.
- Implement BSSN formalism (modern and currently wip).
- Add analytical approximations and tests.
- Improve stability (Spectral BSSN/z4c) and reduce constraints violations.
- Gravitational waves and binary black holes template.

## Contributions
Contributions are welcome! If you wish to contribute, feel free to open a pull request or issue.

## License
GPL-3.0 license. See `LICENSE` for details.

## Contact
For inquiries, reach out via GitHub Issues or email.

---
This project is developed as part of personal advanced numerical relativity research and high-performance computing optimizations. The code and physics computations might have inconsistencies (still at the beginning of the project, and I'm not a physicist).
