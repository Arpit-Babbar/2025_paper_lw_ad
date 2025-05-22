import Pkg
# Activate TenkaiAD environment (TODO - can we include it instead?)
Pkg.activate(joinpath(@__DIR__, "."))
using Tenkai
using TrixiBase

tenkaiad_data_dir = "paper_data"
include(joinpath("plotting", "base_plotting.jl"))
figures_dir = joinpath(@__DIR__, "..", "..", "figures")

mkpath(tenkaiad_data_dir)

isentropic2d = joinpath("isentropic2d")

solver2string = Dict(
    LWEnzymeTower() => "enzyme_tower",
    "lwfr" => "fd",
    MDRKEnzymeTower() => "mdrk_ad",
)

trixi_include_isentropic(; nx, ny, degree, run_number, solver,
                          extra_filename = "", bflux = evaluate) = trixi_include(
    joinpath("examples/run_isentropic.jl"),
    nx = nx, ny = ny,
    saveto = joinpath(
        isentropic2d, solver2string[solver], extra_filename, "$degree", "$nx", "$run_number"),
    degree = degree,
    solver = solver,
    bflux = bflux,
    numerical_flux = Tenkai.EqEuler2D.hllc
)

nx = ny = 50
degree = 1
run_number = 1
solver = LWEnzymeTower()
trixi_include_isentropic(; nx, ny, degree, solver, run_number)

for degree in [1, 2, 3, 4], nx in [20, 40, 80, 160], solver in ["lwfr", LWEnzymeTower()],
    run_number in 1:3
    trixi_include_isentropic(; nx, ny = nx, degree, solver, run_number)
end

# N = 5 for isentropic
for degree in [5], nx in [20, 40, 80, 160], solver in [LWEnzymeTower()], run_number in 1:1
    trixi_include_isentropic(; nx, ny = nx, degree, solver, run_number)
end

for degree in [3], nx in [20, 40, 80, 160], run_number in [1]
    trixi_include_isentropic(; nx, ny = nx, degree, run_number, solver = MDRKEnzymeTower())
end

# Read the data with JSON3.read("output/timer.json")[:inner_timers][Symbol("Cell Residual")][:time_ns] * 1e-9
