using Tenkai
# Submodules
include(joinpath(@__DIR__, "..", "equations", "EqIsentropicEuler1D.jl"))
Eq = EqIsentropicEuler1D
rusanov = Eq.rusanov

using StaticArrays
#------------------------------------------------------------------------------
xmin, xmax = -1.0, 1.0

boundary_condition = (neumann, neumann)
gamma = 1.4
kappa = 1.0

function initial_value_double_rarefaction(x)
    rho = 1000.0
    vel = 3.9
    if x <= 0.0
        v = -vel
    else
        v = vel
    end
    return SVector(rho, rho * v)
end

initial_value = initial_value_double_rarefaction
exact_solution = (x,t) -> initial_value_double_rarefaction(x)

boundary_value = exact_solution # dummy function

degree = 3
solver = LWEnzymeTower()
solution_points = "gl"
correction_function = "radau"
numerical_flux = Eq.rusanov # AD: rusanov_for_ad, non-AD: Eq.rusanov
bound_limit = "yes"
bflux = evaluate
final_time = 0.2

nx = 200
cfl = 0.0
bounds = ([-Inf], [Inf]) # Not used in Euler
tvbM = 0.0
save_iter_interval = 0
save_time_interval = 0.0
animate = true # Factor on save_iter_interval or save_time_interval
compute_error_interval = 1

# blend parameters
indicator_model = "gassner"
debug_blend = false
cfl_safety_factor = 0.95
pure_fv = false

#------------------------------------------------------------------------------
grid_size = nx
domain = [xmin, xmax]
equation = Eq.get_equation(gamma, kappa)
problem = Problem(domain, initial_value, boundary_value, boundary_condition,
                  final_time, exact_solution)
limiter = setup_limiter_blend(blend_type = fo_blend(equation),
                  # indicating_variables = Eq.rho_p_indicator!,
                  indicating_variables = Eq.rho_p_indicator!,
                    # indicating_variables = "conservative",
                  reconstruction_variables = conservative_reconstruction,
                  indicator_model = indicator_model,
                  debug_blend = debug_blend,
                  pure_fv = pure_fv)

# limiter = setup_limiter_tvb(equation; tvbM = 0.0)
scheme = Scheme(solver, degree, solution_points, correction_function,
                numerical_flux, bound_limit, limiter, bflux)
param = Parameters(grid_size, cfl, bounds, save_iter_interval, save_time_interval,
                   compute_error_interval, animate = animate,
                   cfl_safety_factor = 0.5,
                   saveto = "none")
#------------------------------------------------------------------------------
sol = Tenkai.solve(equation, problem, scheme, param);

println(sol["errors"])

return sol;

sol["plot_data"].p_ua