# Simulation_W_JSON

A C++ modular Black Oil Reservoir Simulation project that reads configuration from JSON files.  
Built with C++17, Eigen, nlohmann/json with portable build using CMake. Compatible with MSYS2 (MinGW) or Visual Studio on Windows.

---

## Purpose

This project demonstrates numerical simulation of multiphase fluid flow in reservoirs, using finite difference method.

## Features

- Reads simulation parameters from `resources/config.json` (comes with a config.json example in `resources/`).
- Modular design with headers and source files in `include/` and `src/`.
- Uses **nlohmann/json** for JSON parsing.
- Uses **Eigen** for numerical solvers.
- Portable build with CMake.
- Automatic copying of resource files to the build folder.
- Under `debugging_tools/` there is json_debugger_sat.py, it serves as a novelty visualisation script for the basic simulation variables.

---

## Prerequisites

- **C++17 compiler** (g++ / MSVC)  
- **CMake ≥ 3.16**  
- **MSYS2 (MinGW)** or **Visual Studio** on Windows  

### Install dependencies on MSYS2:

```bash
# Update the system
pacman -Syu
# Close and reopen MSYS2, then update remaining packages
pacman -Su
# Install build tools, compiler, and CMake
pacman -S base-devel mingw-w64-x86_64-toolchain cmake


# Clone the repository
git clone https://github.com/miningcodber/Cpp-Reservoir-Simulator-black-oil-.git
cd Cpp-Reservoir-Simulator-black-oil-

# Create and enter a build folder
mkdir build
cd build

# Configure the project
cmake .. -G "MinGW Makefiles"   # or use your preferred generator

# Build the project
cmake --build .

# Run the program
./Simulation_W_JSON          # uses resources/config.json by default
```
---

#### Overview of core classes and modules

There are 6 core classes used for simulation:

- Grid:
Sets a standard 3D cartesian reservoir, by setting the number of the cells in the x-y-z direction and the size of each cell by setting the cell dimensions in meters (dx, dy, dz). The Grid class contains methods that utilize the row major layout to index specific cells either with coordinates or identification number.

- Rock Properties:
Sets the porosity and permeability factor for each cell (for now both properties are homogeneous in all directions, although heterogeneity can be by applying a Gaussian distribution loop over the private vectors: porosity_ & permeability_).

- Transmissibility Calculator:
Calculates the geometric transmissibility parameter for all directions for each cell, utilizing Grid & Permeability members (it implements harmonic average on permeabilities, calculating the effective permeability for through flow, important on layer stacking).

- Fluid Properties:
The Fluid Properties contain two namespaces, UnitConversion & Correlation. UnitConversion contains inline functions that convert simulation variables from one unit of measurement to another. Correlation also contains inline functions used to calculate PVT variables split upon 4 categories, Gas PVT correlations (BWD 8 constant), Gas Viscosity correlations  (Lee – Gonzalez – Eakin), Oil PVT correlations (Vasquez & Beggs for API<30 & API >30) and Oil Viscosity (Khan et al.). Relative permeability is computed by VERY simplified Corey like curves (no residual  or connate saturations are implemented at this version, nor are taken into consideration in the calculations). Properties are calculated by lambda functions which are used in the final Getters.

- Simulator:
The Simulator utilizes all previous objects made by their constructor respectively. It consists of 4 main methods: initializeState, solveTimestep, computeResidual, assembleJacobian. InitializeState sets the average (homogeneous) water – oil – gas saturation and sets the hydrostatic pressure of the reservoir, calculated by the surface depth. SolveTimestep assembles the newton loop implementing the finite difference method, solves J * dx = -residual, updates variables using newton damping (dx) and checks for convergence. ComputeResidual loops over cells to compute accumulation, implements flag logic for Rs, computes transmissibilities & flux, all to enable positive upwind scheme flow computing. AssembleJacobian assembles the Jacobian matrix, implementing relative epsilon for adaptive perturbation. All methods are tied by the run method.

- Logger:
It’s a class that utilizes nlohmann json to log step#, time, scalar values, pressure, water saturation, gas saturation, and oil saturation. It saves all logged parameters in a .jsonl file under build/.