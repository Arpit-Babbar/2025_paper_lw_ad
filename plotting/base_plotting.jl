using Tenkai: repro_dir
include("$repro_dir/plotting.jl")
rcParams["font.size"] = 16
rcParams["axes.titlesize"] = 16
rcParams["axes.labelsize"] = 16
using DelimitedFiles
using Tenkai.EqLinAdv1D: mult1d
using Tenkai.EqEuler1D: exact_solution_data
import PyPlot as plt
using Plots
using DelimitedFiles
using Plots
using SimpleUnPack: @unpack
using LinearAlgebra
using Printf

markers_array_new() = ["s", "o", "^", "*", "D", ">", "p"]
colors_array_new() = ["red", "royalblue", "green", "m", "c", "y", "k"]
function plot_solns(files::Vector{String}, labels::Vector{String};
                    test_case = nothing, title = "",
                    xlims = nothing, ylims = nothing,
                    exact_data = nothing,
                    xscale = "linear", yscale = "linear",
                    exact_line_width = 1, soln_line_width = 1,
                    linestyles = ["solid" for _ in 1:length(files)],
                    degree = 3,
                    s_ind = 3,
                    s_label = "Pressure",
                    markers = markers_array_new(),
                    colors = colors_array_new(),
                    plt_type = "sol", outdir = nothing)
    @assert length(files) == length(labels)
    fig_density, ax_density = plt.subplots()
    if exact_data === nothing
        exact_data = exact_solution_data(test_case)
    end
    fig_pressure, ax_pressure = plt.subplots()
    if xlims !== nothing
        ax_density.set_xlim(xlims)
        ax_pressure.set_xlim(xlims)
    end
    if ylims !== nothing
        ax_density.set_ylim(ylims)
    end
    ax_density.set_xlabel("x")
    ax_pressure.set_xlabel("x")
    ax_density.set_ylabel("Density")
    ax_pressure.set_ylabel(s_label)
    ax_density.grid(true, linestyle = "--")
    ax_pressure.grid(true, linestyle = "--")

    # Set scales
    ax_density.set_xscale(xscale)
    ax_density.set_yscale(yscale)
    ax_pressure.set_xscale(xscale)
    ax_pressure.set_yscale(yscale)

    if exact_data !== nothing
        @views ax_density.plot(exact_data[:, 1], exact_data[:, 2], label = "Reference",
                               c = "k", linewidth = exact_line_width)
        @views ax_pressure.plot(exact_data[:, 1], exact_data[:, s_ind], label = "Reference",
                                c = "k", linewidth = exact_line_width)
    end
    if plt_type == "avg"
        filename = "avg.txt"
        seriestype = :scatter
    elseif plt_type == "cts_avg"
        filename = "avg.txt"
        seriestype = :line
    else
        @assert plt_type in ("sol", "cts_sol")
        filename = "sol.txt"
        seriestype = :line
    end

    n_plots = length(files)
    for i in 1:n_plots
        # data = datas[i]
        label = labels[i]
        soln_data = readdlm(files[i])
        if plt_type == "avg"
            @views ax_density.plot(soln_data[:, 1], soln_data[:, 2], markers[i],
                                   fillstyle = "none",
                                   c = colors[i], label = label,
                                   linestyle = linestyles[i])
            @views ax_pressure.plot(soln_data[:, 1], soln_data[:, s_ind], markers[i],
                                    fillstyle = "none",
                                    c = colors[i], label = label,
                                    linestyle = linestyles[i])
        elseif plt_type in ("cts_sol", "cts_avg")
            @views ax_density.plot(soln_data[:, 1], soln_data[:, 2], fillstyle = "none",
                                   c = colors[i], label = label,
                                   linewidth = soln_line_width,
                                   linestyle = linestyles[i])
            @views ax_pressure.plot(soln_data[:, 1], soln_data[:, s_ind], fillstyle = "none",
                                    c = colors[i], label = label,
                                    linewidth = soln_line_width,
                                    linestyle = linestyles[i])
        else
            nu = max(2, degree + 1)
            nx = Int(size(soln_data, 1) / nu)
            @views ax_density.plot(soln_data[1:nu, 1], soln_data[1:nu, 2],
                                   fillstyle = "none",
                                   color = colors[i], label = label,
                                   linewidth = soln_line_width,
                                   linestyle = linestyles[i])
            @views ax_pressure.plot(soln_data[1:nu, 1], soln_data[1:nu, s_ind],
                                    fillstyle = "none",
                                    color = colors[i], label = label,
                                    linewidth = soln_line_width,
                                    linestyle = linestyles[i])
            for ix in 2:nx
                i1 = (ix - 1) * nu + 1
                i2 = ix * nu
                @views ax_density.plot(soln_data[i1:i2, 1], soln_data[i1:i2, 2],
                                       fillstyle = "none",
                                       c = colors[i], linewidth = soln_line_width,
                                       linestyle = linestyles[i])
                @views ax_pressure.plot(soln_data[i1:i2, 1], soln_data[i1:i2, s_ind],
                                        fillstyle = "none",
                                        c = colors[i], linewidth = soln_line_width,
                                        linestyle = linestyles[i])
            end
        end
    end
    ax_density.legend()
    ax_pressure.legend()
    ax_density.set_title(title)
    ax_pressure.set_title(title)

    if outdir === nothing
        outdir = joinpath(Tenkai.base_dir, "figures", test_case)
    end
    my_save_fig_python(test_case, fig_density, "density.pdf", fig_dir = outdir)
    my_save_fig_python(test_case, fig_density, "density.png", fig_dir = outdir)
    my_save_fig_python(test_case, fig_pressure, "pressure.pdf", fig_dir = outdir)
    my_save_fig_python(test_case, fig_pressure, "pressure.png", fig_dir = outdir)
    plt.close()
end