# Automatic differentiation for Lax-Wendroff-type discretizations

In order to generate the results from this repository, you need to install `julia`. We recommend using `juliaup`, as detailed in the official website [https://julialang.org](https://julialang.org). The results have been generated using `julia` version 1.10.8, and we recommend installing the same. Once you have installed `julia`, you can clone this repository, enter this directory and start `julia` with the following steps
```shell
git clone https://github.com/Arpit-Babbar/2025_paper_lw_ad.git
cd paper_2025_lw_ad
julia --project=.
```
Then enter the following commands to generate all the results
```julia
julia> import Pkg; Pkg.instantiate() # Does not need to be re-run the next time you enter the REPL
julia> include("generate_all.jl") # Generate all data, postprocess 1D profiles and convergence plots
julia> include("plotting/plot_wct.jl") # See wall clock time performance comparing ALW and AD on screen
```
If you wish to visualize the 2D figure, you need `paraview` and its command line version `pvpython`. Then, in your shell, you can run
```shell
pvpython m2000.py
```
All the figures are now ready and available in the following locations:
1. Double rarefaction test: `paper_figures/isentropic/density.pdf`
2. RHD first Riemann problem: `paper_figures/rhd/density.pdf`
3. Convergence analysis of isentropic vortex test for 2-D compressible Euler's equations: `paper_figures/isentropic_error.pdf`
4. Density profile of Mach 2000 astrophysical jet flow: `paper_figures/m2000.png`.
