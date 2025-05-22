using Tenkai: utils_dir
using Tenkai: cRKSolver
include("$utils_dir/plot_python_solns.jl")
rcParams["font.size"] = 14.65
rcParams["axes.titlesize"] = 16
rcParams["axes.labelsize"] = 16

function add_theo_factors_degrees!(ax, ncells, error, degree, i,
                                   theo_factor_even, theo_factor_odd)
    if degree isa Int64
        d = degree
    else
        d = parse(Int64, degree)
    end
    min_y = minimum(error[1:(end - 1)])
    @show error, min_y
    xaxis = ncells[(end - 1):end]
    slope = d + 1
    if iseven(slope)
        theo_factor = theo_factor_even
    else
        theo_factor = theo_factor_odd
    end
    y0 = theo_factor * min_y
    y = (1.0 ./ xaxis) .^ slope * y0 * xaxis[1]^slope
    markers = ["s", "o", "*", "^"]
    local label_
    if i == 1
        label_ = "\$N+1\$"
    else
        label_ = ""
    end
    # if i == 1
    ax.loglog(xaxis, y,
            #   label = "\$ O(M^{-$(d + 1)})\$",
              label = label_,
              linestyle = "--",
            #   marker = markers[i],
              c = "grey",
              fillstyle = "none")
    # else
    # ax.loglog(xaxis,y, linestyle = "--", c = "grey")
    # end
end

function my_set_ticks!(ax, log_sub, ticks_formatter; dim = 2, base_major = 2.0)
    # Remove scientific notation and set xticks
    # https://stackoverflow.com/a/49306588/3904031

    function anonymous_formatter(y, _)
        if y > 1e-4
            # y_ = parse(Int64, y)
            y_ = Int64(y)
            if dim == 2
                return "\$$y_^2\$"
            else
                return "\$$y_\$"
            end
        else
            return @sprintf "%.4E" y
        end
    end

    formatter = plt.matplotlib.ticker.FuncFormatter(ticks_formatter)
    # (y, _) -> format"{:.4g}".format(int(y)) ) # format"{:.4g}".format(int(y)))
    # https://stackoverflow.com/questions/30887920/how-to-show-minor-tick-labels-on-log-scale-with-plt.matplotlib
    x_major = plt.matplotlib.ticker.LogLocator(base = base_major, subs = (log_sub,),
                                               numticks = 20) # ticklabels at base_major^i*log_sub
    x_minor = plt.matplotlib.ticker.LogLocator(base = 2.0,
                                               subs = LinRange(1.0, 9.0, 9) * 0.1,
                                               numticks = 10)
    #  Used to manipulate tick labels. See help(plt.matplotlib.ticker.LogLocator) for details)
    ax.xaxis.set_major_formatter(anonymous_formatter)
    ax.xaxis.set_minor_formatter(plt.matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_locator(x_major)
    ax.xaxis.set_minor_locator(x_minor)
    ax.tick_params(axis = "both", which = "major")
    ax.tick_params(axis = "both", which = "minor")
end

function plot_python_ndofs_vs_y_degrees(files::Vector{String}, labels::Vector{String},
                                        degrees::Vector{Int}, markers, colors,;
                                        saveto,
                                        theo_factor_even = 0.8, theo_factor_odd = 0.8,
                                        title = nothing, log_sub = "2.5",
                                        error_norm = "l2",
                                        ticks_formatter = format_with_powers,
                                        figsize = (6.4, 4.8),
                                        var_index = 2,
                                        base_major = 2.0,
                                        dim = 1,
                                        plt_type = "dofs")
    # @assert error_type in ["l2","L2"] "Only L2 error for now"
    fig_error, ax_error = plt.subplots(figsize = figsize)
    @assert length(files) == length(markers) == length(colors)
    n_plots = length(files)
    n_labels = length(labels)
    @assert plt_type in ("dofs", "cells")

    # add labels using white colour curves
    for i in 1:n_labels
        data = readdlm(files[i])
        marker = markers[i]
        degree = degrees[i]
        local x
        if plt_type == "dofs"
            x = data[:, 1] * (degree+1)
        else
            x = data[:, 1]
        end
        ax_error.loglog(x, data[:, var_index], marker = marker,
                        # c = "white", # Uncomment this if there are multiple plots with same degree and you want them to be of the same color
                        c = colors[i],
                        mec = "k", fillstyle = "none", label = labels[i])
    end

    for i in 1:n_plots
        data = readdlm(files[i])
        marker = markers[i]
        degree = degrees[i]
        local x
        if plt_type == "dofs"
            x = data[:, 1] * (degree+1)
        else
            x = data[:, 1]
        end
        ax_error.loglog(x, data[:, var_index], marker = marker, c = colors[i],
                        mec = "k", fillstyle = "none")
    end

    for i in eachindex(unique(degrees)) # Assume degrees are not repeated
        data = readdlm(files[i])
        degree = degrees[i]
        local x
        if plt_type == "dofs"
            x = data[:, 1] * (degree+1)
        else
            x = data[:, 1]
        end
        add_theo_factors_degrees!(ax_error, x, data[:,var_index], degree, i,
                                  theo_factor_even, theo_factor_odd)
    end
    if plt_type == "dofs"
        ax_error.set_xlabel("Degrees of freedom")
    else
        ax_error.set_xlabel("Number of elements")
    end
    ax_error.set_ylabel(error_label(error_norm))

    my_set_ticks!(ax_error, log_sub, ticks_formatter; dim = dim, base_major = base_major)

    ax_error.grid(true, linestyle = "--")

    if title !== nothing
        ax_error.set_title(title)
    end
    ax_error.legend()

    fig_error.savefig("$saveto.pdf")
    fig_error.savefig("$saveto.png")

    return fig_error
end

function merge_files(dir)
    ncells = [20, 40, 80, 160]
    errors = []
    for cell in ncells
        file = readdlm(joinpath(dir, "$cell", "1", "error.txt"))
        push!(errors, file[2, 3])
    end
    writedlm(joinpath(dir, "error.txt"), hcat(ncells, errors))
end

file_base = joinpath(@__DIR__, "..", "isentropic2d")

for N in [1, 2, 3, 4, 5]
    merge_files(joinpath(file_base, "enzyme_tower", "$N"))
end

merge_files(joinpath(file_base, "mdrk_ad", "3"))

lw_file(N) = joinpath(file_base, "enzyme_tower", "$N", "error.txt")

mdrk_file = joinpath(file_base, "mdrk_ad", "3", "error.txt")

files = [(lw_file(N) for N in 1:5)..., mdrk_file]
degrees = [1, 2, 3, 4, 5, 3]
labels = ["LW1", "LW2", "LW3", "LW4", "LW5", "MDRK3" ]
markers = ["s", "^", "o", "*",  "D" , ">"]
colors = ["orange", "royalblue", "green", "m", "c" , "y"]
title = "Isentropic vortex"
log_sub = "2.5"
saveto = joinpath("paper_figures", "isentropic_error")
mkpath(saveto)
plot_python_ndofs_vs_y_degrees(files, labels, degrees, markers, colors, title = title,
                               log_sub = log_sub, saveto = saveto,
                               figsize = (6.0, 6.4), plt_type = "cells",
                               dim = 2)
