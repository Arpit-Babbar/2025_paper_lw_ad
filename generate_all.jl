import Pkg
Pkg.activate(joinpath(@__DIR__, "."))
include("generate_1d.jl")
include("generate_isentropic.jl")
include("plotting/plot_convergence_degrees.jl") # Generate convergence plot
include("generate_m2000.jl") # Multiple threads recommended
