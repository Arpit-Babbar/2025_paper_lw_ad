module EqIsentropicEuler1D

using Tenkai.DelimitedFiles
using Tenkai.Plots
using Tenkai.LinearAlgebra
using Tenkai.UnPack
using Tenkai.Printf
using Tenkai.TimerOutputs
using Tenkai.StaticArrays
using Tenkai.Polyester
using Tenkai.LoopVectorization
using Tenkai.JSON3

using Tenkai
using Tenkai.Basis

import Tenkai: admissibility_tolerance

(import Tenkai: flux, prim2con, prim2con!, con2prim, con2prim!,
                eigmatrix,
                limit_slope, zhang_shu_flux_fix,
                apply_tvb_limiter!, apply_bound_limiter!, initialize_plot,
                write_soln!, compute_time_step, post_process_soln,
                flux, update_ghost_values_lwfr!,
                update_ghost_values_fn_blend!,
                correct_variable_bound_limiter!)

(using Tenkai: PlotData, data_dir, get_filename, neumann, minmod,
               get_node_vars,
               set_node_vars!,
               nvariables, eachvariable,
               add_to_node_vars!, subtract_from_node_vars!,
               multiply_add_to_node_vars!,
               update_ghost_values_fn_blend!)

using Tenkai.MuladdMacro

# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

struct IsentropicEuler1D{RealT <: Real} <: AbstractEquations{1, 2}
    gamma::RealT
    kappa::RealT
    nvar::Int
    name::String
    initial_values::Dict{String, Function}
    numfluxes::Dict{String, Function}
end

function flux(x, u, equations::IsentropicEuler1D)
    @unpack gamma = equations
    rho, rho_v = u
    p = pressure(equations, u)
    f1 = rho_v
    f2 = rho_v / rho + p
    return SVector(f1, f2)
end

flux(u, equations::IsentropicEuler1D) = flux(0.0, u, equations)

function density(equations::IsentropicEuler1D, u)
    return u[1]
end

function pressure(equations::IsentropicEuler1D, u)
    rho = u[1]
    p = equations.kappa * (rho^equations.gamma)
    return p
end

function con2prim(equations::IsentropicEuler1D, u)
    rho, rho_v = u
    v = rho_v / rho
    return SVector(rho, v)
end

function con2prim!(eq::IsentropicEuler1D, u_in, u_out)
    u_out .= con2prim(eq, u_in)
end

function prim2con(equations::IsentropicEuler1D, u)
    rho, v = u
    rho_v = rho * v
    return SVector(rho, rho_v)
end

# function prim2con!(eq, u_in, u_out)
#     u_out .= prim2con(eq, u_in)
# end

function max_abs_eigen_value(equations::IsentropicEuler1D, u, dir)
    @unpack gamma = equations
    rho, rho_v = u
    p = pressure(equations, u)
    v = rho_v / rho
    c = sqrt(gamma * p / rho)
    return abs(v) + c
end

function rusanov(x, ual, uar, Fl, Fr, Ul, Ur, eq::IsentropicEuler1D, dir)
    @unpack gamma = eq
    lamba_ll = max_abs_eigen_value(eq, ual, dir)
    lamba_rr = max_abs_eigen_value(eq, uar, dir)
    λ = max(lamba_ll, lamba_rr) # local wave speed
    f1 = 0.5 * (Fl[1] + Fr[1]) - 0.5 * λ * (Ur[1] - Ul[1])
    f2 = 0.5 * (Fl[2] + Fr[2]) - 0.5 * λ * (Ur[2] - Ul[2])
    return SVector(f1, f2)
end

rusanov_for_ad(x, ual, uar, Fl, Fr, Ul, Ur, eq::IsentropicEuler1D, dir) =
    rusanov(x, ual, uar, Fl, Fr, Ul, Ur, eq, dir)

function compute_time_step(eq::IsentropicEuler1D, problem, grid, aux, op, cfl, u1, ua)
    @unpack source_terms = problem
    nx = grid.size
    dx = grid.dx
    den = 0.0
    for i in 1:nx
        u = get_node_vars(ua, eq, i)
        smax = max_abs_eigen_value(eq, u, 1)
        den = max(den, smax / dx[i])
    end
    dt = cfl / den
    return dt
end

function correct_variable_bound_limiter!(variable::typeof(density), eq::AbstractEquations{1},
                                         grid, op, ua, u1)
    @unpack Vl, Vr = op
    nx = grid.size
    nd = op.degree + 1
    eps = 1e-10
    for element in 1:nx
        var_ll = var_rr = 0.0
        var_min = 1e20
        for i in Base.OneTo(nd)
            u_node = get_node_vars(u1, eq, i, element)
            var = variable(eq, u_node)
            var_ll += var * Vl[i]
            var_rr += var * Vr[i]
            var_min = min(var_min, var)
        end
        var_min = min(var_min, var_ll, var_rr)
        ua_ = get_node_vars(ua, eq, element)
        var_avg = variable(eq, ua_)
        @assert var_avg>0.0 "Failed at element $element", var_avg
        eps_ = min(eps, 0.1 * var_avg)
        ratio = abs(eps_ - var_avg) / (abs(var_min - var_avg) + 1e-13)
        theta = min(ratio, 1.0) # theta for preserving positivity of density
        if theta < 1.0
            for i in 1:nd
                u_node = get_node_vars(u1, eq, i, element)
                multiply_add_set_node_vars!(u1,
                                            theta, u_node,
                                            1 - theta, ua_,
                                            eq, i, element)
            end
        end
    end
end

function Tenkai.apply_bound_limiter!(eq::IsentropicEuler1D, grid, scheme, param, op, ua,
                                     u1, aux)
    if scheme.bound_limit == "no"
        return nothing
    end
    correct_variable_bound_limiter!(density, eq, grid, op, ua, u1)
end

@inbounds @inline function rho_p_indicator!(un, eq::IsentropicEuler1D)
    for ix in 1:size(un, 2) # loop over dofs and faces
        u_node = get_node_vars(un, eq, ix)
        p = pressure(eq, u_node)
        un[1, ix] *= p # ρ * p
    end
    n_ind_var = 1
    return n_ind_var
end

function Tenkai.zhang_shu_flux_fix(eq::IsentropicEuler1D,
                                   uprev,    # Solution at previous time level
                                   ulow,     # low order update
                                   Fn,       # Blended flux candidate
                                   fn_inner, # Inner part of flux
                                   fn,       # low order flux
                                   c)
    uhigh = uprev - c * (Fn - fn_inner) # First candidate for high order update
    rho_low, rho_high = density(eq, ulow), density(eq, uhigh)
    eps = 0.1 * rho_low
    ratio = abs(eps - rho_low) / (abs(rho_high - rho_low) + 1e-13)
    theta = min(ratio, 1.0)
    if theta < 1.0
        Fn = theta * Fn + (1.0 - theta) * fn # Second candidate for flux
    end
end

function Tenkai.limit_slope(eq::IsentropicEuler1D, slope, ufl, u_star_ll, ufr, u_star_rr,
                            ue, xl, xr, el_x = nothing, el_y = nothing)

    # The MUSCL-Hancock scheme is guaranteed to be admissibility preserving if
    # slope is chosen so that
    # u_star_l = ue + 2.0*slope*xl, u_star_r = ue+2.0*slope*xr are admissible
    # ue is already admissible and we know we can find sequences of thetas
    # to make theta*u_star_l+(1-theta)*ue is admissible.
    # This is equivalent to replacing u_star_l by
    # u_star_l = ue + 2.0*theta*s*xl.
    # Thus, we simply have to update the slope by multiplying by theta.

    slope, u_star_ll, u_star_rr = limit_variable_slope(eq, density, slope,
                                                       u_star_ll, u_star_rr, ue, xl, xr)

    ufl = ue + slope * xl
    ufr = ue + slope * xr

    return ufl, ufr, slope
end

varnames(::Type{IsentropicEuler1D}) = ["ρ", "ρv"]
varnames(::Type{IsentropicEuler1D}, i::Int) = varnames(IsentropicEuler1D)[i]

function eigmatrix(eq::IsentropicEuler1D, u)
    Id = SMatrix{nvariables(eq), nvariables(eq)}(1.0, 0.0,
                                                 0.0, 1.0)

    return Id, Id
end

function Tenkai.apply_tvb_limiter!(eq::IsentropicEuler1D, problem, scheme, grid, param, op, ua,
                                   u1, aux)
    @timeit aux.timer "TVB limiter" begin
    #! format: noindent
    nx = grid.size
    @unpack xg, wg, Vl, Vr = op
    @unpack limiter = scheme
    @unpack tvbM, cache = limiter
    left_bc, right_bc = problem.boundary_condition
    nd = length(wg)
    nvar = nvariables(eq)
    # face values
    (uimh, uiph, Δul, Δur, Δual, Δuar, char_Δul, char_Δur, char_Δual, char_Δuar,
    dulm, durm, du) = cache

    # Loop over cells
    beta = limiter.beta
    beta_ = beta / 2.0
    for cell in 1:nx
        ual, ua_, uar = (get_node_vars(ua, eq, cell - 1),
                         get_node_vars(ua, eq, cell),
                         get_node_vars(ua, eq, cell + 1))

        # Needed for characteristic limiting
        R, L = eigmatrix(eq, ua_)
        fill!(uimh, zero(eltype(uimh)))
        fill!(uiph, zero(eltype(uiph)))
        Mdx2 = tvbM * grid.dx[cell]^2
        if left_bc == neumann && right_bc == neumann && (cell == 1 || cell == nx)
            Mdx2 = 0.0 # Force TVD on boundary for Shu-Osher
        end
        # end # timer
        for ii in 1:nd
            u_ = get_node_vars(u1, eq, ii, cell)
            multiply_add_to_node_vars!(uimh, Vl[ii], u_, eq, 1)
            multiply_add_to_node_vars!(uiph, Vr[ii], u_, eq, 1)
        end
        # Get views of needed cell averages
        # slopes b/w centres and faces

        uimh_ = get_node_vars(uimh, eq, 1)
        uiph_ = get_node_vars(uiph, eq, 1)

        # We will set
        # Δul[n] = ua_[n] - uimh[n]
        # Δur[n] = uiph[n] - ua_[n]
        # Δual[n] = ua_[n] - ual[n]
        # Δuar[n] = uar[n] - ua_[n]

        set_node_vars!(Δul, ua_, eq, 1)
        set_node_vars!(Δur, uiph_, eq, 1)
        set_node_vars!(Δual, ua_, eq, 1)
        set_node_vars!(Δuar, uar, eq, 1)

        subtract_from_node_vars!(Δul, uimh_, eq)
        subtract_from_node_vars!(Δur, ua_, eq)
        subtract_from_node_vars!(Δual, ual, eq)
        subtract_from_node_vars!(Δuar, ua_, eq)

        Δul_ = get_node_vars(Δul, eq, 1)
        Δur_ = get_node_vars(Δur, eq, 1)
        Δual_ = get_node_vars(Δual, eq, 1)
        Δuar_ = get_node_vars(Δuar, eq, 1)

        # Uncomment this part for characteristic limiting
        # mul!(char_Δul, L, Δul_)   # char_Δul = L*Δul
        # mul!(char_Δur, L, Δur_)   # char_Δur = L*Δur
        # mul!(char_Δual, L, Δual_) # char_Δual = L*Δual
        # mul!(char_Δuar, L, Δuar_) # char_Δuar = L*Δuar

        # Use primitive variables
        # set_node_vars!(char_Δul, con2prim(eq, Δul_), eq, 1)
        # set_node_vars!(char_Δur, con2prim(eq, Δur_), eq, 1)
        # set_node_vars!(char_Δual, con2prim(eq, Δual_), eq, 1)
        # set_node_vars!(char_Δuar, con2prim(eq, Δuar_), eq, 1)

        # Keep conservative variables
        set_node_vars!(char_Δul, Δul_, eq, 1)
        set_node_vars!(char_Δur, Δur_, eq, 1)
        set_node_vars!(char_Δual, Δual_, eq, 1)
        set_node_vars!(char_Δuar, Δuar_, eq, 1)


        char_Δul_ = get_node_vars(char_Δul, eq, 1)
        char_Δur_ = get_node_vars(char_Δur, eq, 1)
        char_Δual_ = get_node_vars(char_Δual, eq, 1)
        char_Δuar_ = get_node_vars(char_Δuar, eq, 1)
        for n in eachvariable(eq)
            dulm[n] = minmod(char_Δul_[n], beta_ * char_Δual_[n], beta_ * char_Δuar_[n], Mdx2)
            durm[n] = minmod(char_Δur_[n], beta_ * char_Δual_[n], beta_ * char_Δuar_[n], Mdx2)
        end

        # limit if jumps are detected
        dulm_ = get_node_vars(dulm, eq, 1)
        durm_ = get_node_vars(durm, eq, 1)
        jump_l = jump_r = 0.0
        for n in 1:nvar
            jump_l += abs(char_Δul_[n] - dulm_[n])
            jump_r += abs(char_Δur_[n] - durm_[n])
        end
        jump_l /= nvar
        jump_r /= nvar

        if jump_l > 1e-06 || jump_r > 1e-06
            add_to_node_vars!(durm, dulm_, eq, 1) # durm = durm + dulm
            # We want durm = 0.5 * (dul + dur), we adjust 0.5 later

            # Uncomment this for characteristic variables
            # mul!(du, R, durm)            # du = R * (dulm+durm)

            # Conservative / primitive variables
            # durm_ = prim2con(eq, get_node_vars(durm, eq, 1)) # Keep for primitive
            durm_ = get_node_vars(durm, eq, 1) # Keep for conservative
            set_node_vars!(du, durm_, eq, 1)
            for ii in Base.OneTo(nd)
                du_ = get_node_vars(du, eq, 1)
                set_node_vars!(u1, ua_ + (xg[ii] - 0.5) * du_, # 2.0 adjusted with 0.5 above
                               eq, ii,
                               cell)
            end
        end
    end
    return nothing
    end # timer
end

#-------------------------------------------------------------------------------
# Plotting functions
#-------------------------------------------------------------------------------

function Tenkai.initialize_plot(eq::IsentropicEuler1D, op, grid, problem, scheme, timer, u1,
                                ua)
    @timeit timer "Write solution" begin
    #! format: noindent
    @timeit timer "Initialize write solution" begin
    #! format: noindent
    # Clear and re-create output directory
    rm("output", force = true, recursive = true)
    mkdir("output")

    xc = grid.xc
    nx = grid.size
    @unpack xg = op
    nd = op.degree + 1
    nu = max(nd, 2)
    xu = LinRange(0.0, 1.0, nu)
    Vu = Vandermonde_lag(xg, xu)
    xf = grid.xf
    nvar = eq.nvar
    # Create plot objects to be later collected as subplots

    # Creating a subplot for title
    p_title = plot(title = "Cell averages plot, $nx cells, t = 0.0",
                   grid = false, showaxis = false, bottom_margin = 0Plots.px)
    # Initialize subplots for density, velocity and pressure
    p_ua, p_u1 = [plot() for _ in 1:nvar], [plot() for _ in 1:nvar]
    labels = ["Density", "Velocity", "Pressure"]
    y = zeros(nx) # put dummy to fix plotly bug with OffsetArrays
    for n in 1:nvar
        @views plot!(p_ua[n], xc, y, label = "Approximate",
                     linestyle = :dot, seriestype = :scatter,
                     color = :blue, markerstrokestyle = :dot,
                     markershape = :circle, markersize = 2,
                     markerstrokealpha = 0)
        xlabel!(p_ua[n], "x")
        ylabel!(p_ua[n], labels[n])
    end
    l_super = @layout[a{0.01h}; b c] # Selecting layout for p_title being title
    p_ua = plot(p_title, p_ua[1], p_ua[2], layout = l_super,
                size = (1000, 500)) # Make subplots

    # Set up p_u1 to contain polynomial approximation as a different curve
    # for each cell
    x = LinRange(xf[1], xf[2], nu)
    up1 = zeros(nvar, nd)
    u = zeros(nu)
    for ii in 1:nd
        @views con2prim!(eq, u1[:, ii, 1], up1[:, ii]) # store prim form in up1
    end

    for n in 1:nvar
        u = @views Vu * up1[n, :]
        plot!(p_u1[n], x, u, color = :red, legend = false)
        xlabel!(p_u1[n], "x")
        ylabel!(p_u1[n], labels[n])
    end
    for i in 2:nx
        for ii in 1:nd
            @views con2prim!(eq, u1[:, ii, i], up1[:, ii]) # store prim form in up1
        end
        x = LinRange(xf[i], xf[i + 1], nu)
        for n in 1:nvar
            u = @views Vu * up1[n, :]
            plot!(p_u1[n], x, u, color = :red, label = nothing, legend = false)
        end
    end

    l = @layout[a{0.01h}; b c] # Selecting layout for p_title being title
    p_u1 = plot(p_title, p_u1[1], p_u1[2], layout = l,
                size = (1200, 500)) # Make subplots

    anim_ua, anim_u1 = Animation(), Animation() # Initialize animation objects
    plot_data = PlotData(p_ua, anim_ua, p_u1, anim_u1)
    return plot_data
    end # timer
    end # timer
end

function Tenkai.write_soln!(base_name, fcount, iter, time, dt, eq::IsentropicEuler1D, grid,
                            problem, param, op, ua, u1, aux, ndigits = 3)
    @timeit aux.timer "Write solution" begin
    #! format: noindent
    @unpack plot_data = aux
    avg_filename = get_filename("output/avg", ndigits, fcount)
    @unpack p_ua, p_u1, anim_ua, anim_u1 = plot_data
    @unpack final_time = problem
    xc = grid.xc
    nx = grid.size
    @unpack xg = op
    nd = op.degree + 1
    nu = max(nd, 2)
    xu = LinRange(0.0, 1.0, nu)
    Vu = Vandermonde_lag(xg, xu)
    nvar = eq.nvar
    @unpack save_time_interval, save_iter_interval, animate = param
    avg_file = open("$avg_filename.txt", "w")
    up_ = zeros(nvar)
    ylims = [[Inf, -Inf] for _ in 1:nvar] # set ylims for plots of all variables
    for i in 1:nx
        @views con2prim!(eq, ua[:, i], up_) # store primitve form in up_
        @printf(avg_file, "%e %e %e\n", xc[i], up_[1], up_[2])
        # TOTHINK - Check efficiency of printf
        for n in 1:(eq.nvar)
            p_ua[n + 1][1][:y][i] = @views up_[n]    # Update y-series
            ylims[n][1] = min(ylims[n][1], up_[n]) # Compute ymin
            ylims[n][2] = max(ylims[n][2], up_[n]) # Compute ymax
        end
    end
    close(avg_file)
    for n in 1:nvar # set ymin, ymax for ua, u1 plots
        ylims!(p_ua[n + 1], (ylims[n][1] - 0.1, ylims[n][2] + 0.1))
        ylims!(p_u1[n + 1], (ylims[n][1] - 0.1, ylims[n][2] + 0.1))
    end
    t = round(time; digits = 3)
    title!(p_ua[1], "Cell averages plot, $nx cells, t = $t")
    sol_filename = get_filename("output/sol", ndigits, fcount)
    sol_file = open(sol_filename * ".txt", "w")
    up1 = zeros(nvar, nd)

    u = zeros(nvar, nu)
    x = zeros(nu)
    for i in 1:nx
        for ii in 1:nd
            @views con2prim!(eq, u1[:, ii, i], up1[:, ii]) # store prim form in up1
        end
        @. x = grid.xf[i] + grid.dx[i] * xu
        @views mul!(u, up1, Vu')
        for n in 1:nvar
            p_u1[n + 1][i][:y] = u[n, :]
        end
        for ii in 1:nu
            @printf(sol_file, "%e %e %e\n", x[ii], u[1, ii], u[2, ii])
        end
    end
    close(sol_file)
    title!(p_u1[1], "Numerical Solution, $nx cells, t = $t")
    println("Wrote $sol_filename.txt, $avg_filename.txt")
    if problem.final_time - time < 1e-10
        cp("$avg_filename.txt", "./output/avg.txt", force = true)
        cp("$sol_filename.txt", "./output/sol.txt", force = true)
        println("Wrote final solution to avg.txt, sol.txt.")
    end
    if animate == true
        if abs(time - final_time) < 1.0e-10
            frame(anim_ua, p_ua)
            frame(anim_u1, p_u1)
        end
        if save_iter_interval > 0
            animate_iter_interval = save_iter_interval
            if mod(iter, animate_iter_interval) == 0
                frame(anim_ua, p_ua)
                frame(anim_u1, p_u1)
            end
        elseif save_time_interval > 0
            animate_time_interval = save_time_interval
            k1, k2 = ceil(time / animate_time_interval),
                     floor(time / animate_time_interval)
            if (abs(time - k1 * animate_time_interval) < 1e-10 ||
                abs(time - k2 * animate_time_interval) < 1e-10)
                frame(anim_ua, p_ua)
                frame(anim_u1, p_u1)
            end
        end
    end
    fcount += 1
    return fcount
    end # timer
end

function exact_solution_data(test_case)
    @warn "Exact solution does not set!"
    return nothing
end

function Tenkai.post_process_soln(eq::IsentropicEuler1D, aux, problem, param, scheme)
    @unpack timer, error_file = aux
    @timeit timer "Write solution" begin
    #! format: noindent
    println("Post processing solution")
    nvar = eq.nvar
    @unpack plot_data = aux
    @unpack p_ua, p_u1, anim_ua, anim_u1 = plot_data
    @unpack animate, saveto = param
    initial_values = eq.initial_values
    if problem.initial_value in values(initial_values) # Using ready made tests
        initial_value_string, = [a
                                 for (a, b) in initial_values
                                 if
                                 b == problem.initial_value]
        exact_data = exact_solution_data(initial_value_string)

        for n in 1:nvar
            @views plot!(p_ua[n + 1], exact_data[:, 1], exact_data[:, n + 1],
                         label = "Exact",
                         color = :black)
            @views plot!(p_u1[n + 1], exact_data[:, 1], exact_data[:, n + 1],
                         label = "Exact",
                         color = :black, legend = true)
            ymin = min(minimum(p_ua[n + 1][1][:y]), minimum(exact_data[:, n + 1]))
            ymax = max(maximum(p_ua[n + 1][1][:y]), maximum(exact_data[:, n + 1]))
            ylims!(p_ua[n + 1], (ymin - 0.1, ymax + 0.1))
            ylims!(p_u1[n + 1], (ymin - 0.1, ymax + 0.1))
        end
    end
    savefig(p_ua, "output/avg.png")
    savefig(p_u1, "output/sol.png")
    savefig(p_ua, "output/avg.html")
    savefig(p_u1, "output/sol.html")
    if animate == true
        gif(anim_ua, "output/avg.mp4", fps = 5)
        gif(anim_u1, "output/sol.mp4", fps = 5)
    end
    println("Wrote avg, sol in gif,html,png format to output directory.")
    plot(p_ua)
    plot(p_u1)

    close(error_file)
    if saveto != "none"
        if saveto[end] == "/"
            saveto = saveto[1:(end - 1)]
        end
        mkpath(saveto)
        for file in readdir("./output")
            cp("./output/$file", "$saveto/$file", force = true)
        end
        cp("./error.txt", "$saveto/error.txt", force = true)
        println("Saved output files to $saveto")
    end
    end # timer

    # Print timer data on screen
    print_timer(aux.timer, sortby = :firstexec)
    print("\n")
    show(aux.timer)
    print("\n")
    println("Time outside write_soln = "
            *
            "$(( TimerOutputs.tottime(timer)
                - TimerOutputs.time(timer["Write solution"]) ) * 1e-9)s")
    println("─────────────────────────────────────────────────────────────────────────────────────────")
    timer_file = open("./output/timer.json", "w")
    JSON3.write(timer_file, TimerOutputs.todict(timer))
    close(timer_file)
    return nothing
end

function get_equation(gamma, kappa)
    name = "Isentropic Euler 1D"
    nvar = 2
    numfluxes = Dict{String, Function}()
    initial_values = Dict{String, Function}()
    return IsentropicEuler1D(gamma, kappa, nvar, name, initial_values, numfluxes)
end

function update_ghost_values_lwfr!(problem, scheme, eq::IsentropicEuler1D, grid, aux, op, cache,
                                   t, dt, scaling_factor = 1)
    @timeit aux.timer "Update ghost values" begin
    #! format: noindent
    @unpack Fb, Ub = cache
    update_ghost_values_periodic!(eq, problem, Fb, Ub)

    if problem.periodic_x
        return nothing
    end

    nx = grid.size
    @unpack degree, xg, wg = op
    nd = degree + 1
    dx, xf = grid.dx, grid.xf
    nvar = nvariables(eq)
    @unpack boundary_value, boundary_condition = problem
    left, right = boundary_condition
    refresh!(u) = fill!(u, 0.0)

    ub, fb = zeros(nvar), zeros(nvar)

    # For Dirichlet bc, use upwind flux at faces by assigning both physical
    # and ghost cells through the bc.
    dt_scaled = scaling_factor * dt
    wg_scaled = scaling_factor * wg # It is a static array so doesn't cause allocations
    if left == dirichlet
        x = xf[1]
        for l in 1:nd
            tq = t + xg[l] * dt_scaled
            ubvalue = boundary_value(x, tq)
            fbvalue = flux(x, ubvalue, eq)
            for n in 1:nvar
                ub[n] += ubvalue[n] * wg_scaled[l]
                fb[n] += fbvalue[n] * wg_scaled[l]
            end
        end
        for n in 1:nvar
            Ub[n, 1, 1] = Ub[n, 2, 0] = ub[n]
            Fb[n, 1, 1] = Fb[n, 2, 0] = fb[n]
        end
    elseif left == neumann
        for n in 1:nvar
            Ub[n, 2, 0] = Ub[n, 1, 1]
            Fb[n, 2, 0] = Fb[n, 1, 1]
        end
    elseif left == reflect
        # velocity reflected back in opposite direction and density is same
        for n in 1:nvar
            Ub[n, 2, 0] = Ub[n, 1, 1]
            Fb[n, 2, 0] = Fb[n, 1, 1]
        end
        Ub[2, 2, 0] = -Ub[2, 2, 0] # velocity reflected back
        Fb[1, 2, 0] = -Fb[1, 2, 0] # vel multiple term
    else
        println("Incorrect bc specified at left.")
        @assert false
    end

    refresh!.((ub, fb))
    if right == dirichlet
        x = xf[nx + 1]
        for l in 1:nd
            tq = t + xg[l] * dt_scaled
            ubvalue = boundary_value(x, tq)
            fbvalue = flux(x, ub, eq)
            for n in 1:nvar
                ub[n] += ubvalue[n] * wg_scaled[l]
                fb[n] += fbvalue[n] * wg_scaled[l]
            end
        end
        for n in 1:nvar
            Ub[n, 2, nx] = Ub[n, 1, nx + 1] = ub[n]
            Fb[n, 2, nx] = Fb[n, 1, nx + 1] = fb[n]
        end
    elseif right == neumann
        for n in 1:nvar
            Ub[n, 1, nx + 1] = Ub[n, 2, nx]
            Fb[n, 1, nx + 1] = Fb[n, 2, nx]
        end
    elseif right == reflect
        # velocity reflected back in opposite direction and density is same
        for n in 1:nvar
            Ub[n, 1, nx + 1] = Ub[n, 2, nx]
            Fb[n, 1, nx + 1] = Fb[n, 2, nx]
        end
        Ub[2, 1, nx + 1] = -Ub[2, 1, nx + 1] # velocity reflected back
        Fb[1, 1, nx + 1] = -Fb[1, 1, nx + 1] # vel multiple term

    else
        println("Incorrect bc specified at right.")
        @assert false
    end

    if scheme.limiter.name == "blend"
        update_ghost_values_fn_blend!(eq, problem, grid, aux)
    end
    end # timer
    return nothing
end

end # @muladd

end # module EqEuler1D