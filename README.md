# csc2521 project

The code implements the paper *Vertex Block Descent* by Chen et al. Specifically, the code implements "Algorithm 1: VBD simulation for one time step" in the paper, without discrete and continuous collision detection and without parallelization using vertex graph coloring. The code demonstrates the algorithm with a cantilever beam fixed at one end.

The code relies on Eigen for numerical linear algebra operations, libigl for reading mesh files, and polyscope for visualization.

## How to build

1. Clone the repository and all its submodules recursively.
2. Create a subdirectory `build` in the repository.
3. Run `cd build` and `cmake .. -DCMAKE_BUILD_TYPE=Release` to build in Release mode for runtime speed.
4. Run `make -j <job>` to compile the source code.
5. Run the compiled executable. The mesh file is provided in the repository and is retrieved at the relative path `../data/bar_test.mesh`. The user needs to ensure that the relative path is valid when running the executable.

## Platform

The code has been successfully compiled on Windows 11 with Visual Studio 2022, and on Ubuntu 22.02 with g++.
