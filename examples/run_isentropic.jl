using Tenkai
Eq = Tenkai.EqEuler2D
#------------------------------------------------------------------------------
xmin, xmax = -10.0, 10.0
ymin, ymax = -10.0, 10.0

boundary_value = Eq.zero_bv # dummy
boundary_condition = (periodic, periodic, periodic, periodic)
γ = 1.4

initial_value = Eq.isentropic_iv
exact_solution = Eq.isentropic_exact

degree = 5
solver = LWEnzymeTower()
solution_points = "gl"
correction_function = "radau"
numerical_flux = Eq.rusanov # AD: rusanov_for_ad, non-AD: Eq.rusanov
bound_limit = "no"
bflux = evaluate
final_time = 1.0

nx = 30
ny = 30
cfl = 0.0
bounds = ([-Inf], [Inf]) # Not used in Euler
tvbM = 0.0
save_iter_interval = 0
save_time_interval = 0.0 # final_time / 5.0
animate = true # Factor on save_iter_interval or save_time_interval
compute_error_interval = 0

cfl_safety_factor = 0.8

#------------------------------------------------------------------------------
grid_size = [nx, ny]
domain = [xmin, xmax, ymin, ymax]
equation = Eq.get_equation(γ)
problem = Problem(domain, initial_value, boundary_value, boundary_condition,
                  final_time, exact_solution)
limiter = setup_limiter_none()
scheme = Scheme(solver, degree, solution_points, correction_function,
                numerical_flux, bound_limit, limiter, bflux,
                2)
param = Parameters(grid_size, cfl, bounds, save_iter_interval,
                   save_time_interval, compute_error_interval,
                   animate = animate,
                   saveto = "none", cfl_safety_factor = cfl_safety_factor)
#------------------------------------------------------------------------------
sol = Tenkai.solve(equation, problem, scheme, param);

println(sol["errors"])

return sol;

# using TimerOutputs: todict
# todict(sol["aux"].timer)["inner_timers"]["Cell Residual"]["time_ns"] / 1e9
# Use this to extract only cell residual cost
# TimerOutputs.todict(timer)

# 12.1 sec

# 16.8 sec
# 14.7
