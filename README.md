# Thesis project

## Setup
- make sure you have Julia installed: ([https://julialang.org/downloads/](https://julialang.org/downloads/))
- next install needed packages:
    - open julia in terminal
    ```bash
    julia
    ```
    - install packages
    ```julia
    using Pkg
    Pkg.add("Plots")
    Pkg.add("LinearAlgebra")
    ```

## Run
in the root directory of the project run:
```bash
julia main.jl
```

## Parameters
All parameters are set in the `const.jl` file. The parameters are:
- `N`: number of gaussians
- `M`: mass of particles
- `A`: parameter of the gaussian
- `L`: length of the box
