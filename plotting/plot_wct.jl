using DelimitedFiles
using PrettyTables
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Tenkai: JSON3

isentropic2d = "isentropic2d"

function file2wct(file)
    data = JSON3.read(file)
    total_time = data[:total_time_ns]
    data_inner = data[:inner_timers]
    extra_time = (data_inner[Symbol("Write solution")][:time_ns]
                  + data_inner[Symbol("Limiter")][:time_ns])
    return (total_time - extra_time) * 1e-9
end

base_file(solver, degree, nx) = joinpath(isentropic2d, solver, "$degree", "$nx")

function trunc_float(x)
    str = "$(x)"
    if length(str) <= 4
        return str
    end
    x < 1e-2 ? str[1:5] : str[1:4]
end

function get_wct(;nx, degree, solver)
    base = base_file(solver, degree, nx)
    files = [joinpath(base, "$i", "timer.json") for i in 1:3]
    wcts = file2wct.(files)
    min_val = minimum(wcts)
    if min_val < 1e-2
        "$(minimum(wcts))"[1:5]
    else
        "$(minimum(wcts))"[1:4]
    end
end

function get_wct_vec(; nx_array, degree, solver)
    [get_wct(nx = nx, degree = degree, solver = solver) for nx in nx_array]
end

A = Vector{Any}()
header = [" "]
finest_lw = []
finest_enzyme = []
finest_lwtde = []
for degree in 1:4
    nx_array = [20, 40, 80, 160]
    # efficient_arr = get_wct_vec(; nx_array, degree = "$degree", solver = "efficient")
    fd_arr = get_wct_vec(; nx_array, degree = degree, solver = "fd")
    enzyme_arr = get_wct_vec(; nx_array, degree = degree, solver = "enzyme_tower")
    lwtde_arr = get_wct_vec(; nx_array, degree = degree, solver = "lwtde")
    push!(finest_lw, fd_arr[end])
    push!(finest_enzyme, enzyme_arr[end])
    push!(finest_lwtde, lwtde_arr[end])
    push!(A, fd_arr, enzyme_arr, lwtde_arr)
    push!(header, ("FD", "ADP", "ADE")...)
end

str2float(x) = parse(Float64, x)

fd_arr_float = str2float.(finest_lw)
enzyme_arr_float = str2float.(finest_enzyme)
lwtde_arr_float = str2float.(finest_lwtde)

ratio_fd, ratio_enzyme, ratio_lwtde = (trunc_float.(fd_arr_float ./ fd_arr_float),
          trunc_float.(fd_arr_float ./ enzyme_arr_float),
          trunc_float.(fd_arr_float ./ lwtde_arr_float))

push!(A, fd_arr_float, enzyme_arr_float, lwtde_arr_float)

all_ratios = [ratio_fd..., ratio_enzyme..., ratio_lwtde...]

push!(header, ("FD", "ADP", "ADE")...)

get_wct_vec(; nx_array = [20, 40, 80, 160], degree = 4, solver = "lwtde")

A = hcat(["20", "40", "80", "160"], A...)

pretty_table(A, header = header, backend = Val(:latex))

println("\nSpeedup ratios (FD over others) for finest grid (N=4, 160x160):")

for i in eachindex(all_ratios)
    print(" & ")
    print(all_ratios[i])
    # if i < length(all_ratios)

    # else
    #     println()
    # end
end
