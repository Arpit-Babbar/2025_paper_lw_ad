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
for degree in 1:4
    nx_array = [20, 40, 80, 160]
    # efficient_arr = get_wct_vec(; nx_array, degree = "$degree", solver = "efficient")
    fd_arr = get_wct_vec(; nx_array, degree = degree, solver = "fd")
    enzyme_arr = get_wct_vec(; nx_array, degree = degree, solver = "enzyme_tower")
    push!(A, fd_arr, enzyme_arr)
    push!(header, ("ALW", "AD")...)
end

A = hcat(["20", "40", "80", "160"], A...)

pretty_table(A, header = header, backend = Val(:latex))
