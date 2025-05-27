import Pkg
Pkg.activate(joinpath(@__DIR__, "."))
using Tenkai
using TrixiBase

tenkaiad_data_dir = joinpath("paper_data")
include(joinpath("plotting", "base_plotting.jl"))
figures_dir = joinpath(@__DIR__, "paper_figures")

mkpath(tenkaiad_data_dir)

isentropic_dir = joinpath(tenkaiad_data_dir, "isentropic_euler_double_rarefaction")
isentropic_ref_dir = joinpath(tenkaiad_data_dir, "isentropic_euler_double_rarefaction_ref")
trixi_include("examples/run_isentropic_euler_double_rarefaction.jl", nx = 200, saveto = isentropic_dir)

trixi_include(joinpath("examples/run_isentropic_euler_double_rarefaction.jl"),
                       nx = 50000, saveto = isentropic_ref_dir,
                       degree = 0,
                       solver = "lwfr",
                       bound_limit = "no",
                       limiter = Tenkai.setup_limiter_none(),
                       compute_error_interval = 0)

files = [joinpath(isentropic_ref_dir, "sol.txt"), joinpath(isentropic_dir, "sol.txt")]
labels = ["Reference", "LW3"]

plot_solns(files, labels, outdir = joinpath(figures_dir, "isentropic"),
           linestyles = ["solid", "dotted"],
           colors = ["black", "red"],
           exact_line_width = 2.5, soln_line_width = 2.5,
           plt_type = "cts_avg")

Tenkai.Enzyme.API.strictAliasing!(false)

rhd_dir = "RHD1D_RyuRP_P1"
rhd_ref_dir = "RHD1D_RyuRP_P1_ref"
trixi_include(joinpath("examples/run_RHD1D_RyuRP_P1.jl"), nx = 300, saveto = rhd_dir)

# This one takes time to run!
trixi_include(joinpath("examples/run_RHD1D_RyuRP_P1.jl"), nx = 20000,
              degree = 0,
              limiter = Tenkai.setup_limiter_none(),
              solver = "lwfr",
              bound_limit = "no",
              saveto = rhd_ref_dir,
              animate = false)

files = [joinpath(rhd_ref_dir, "sol.txt"), joinpath(rhd_dir, "sol.txt")]
labels = ["Reference", "LW3"]
plot_solns(files, labels, outdir = joinpath(figures_dir, "rhd"),
            linestyles = ["solid", "dotted"],
            # plt_type = "cts_avg",
            exact_line_width = 2.5, soln_line_width = 2.5,
            colors = ["black", "red"])
