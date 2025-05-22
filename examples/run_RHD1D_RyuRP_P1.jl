using Tenkai
using Tenkai.Plots
include(joinpath(@__DIR__, "..", "equations", "EqRHD1D_geos_1v_comp.jl"))

# Submodules
Eq = EqRHD1D_geos_1v_comp
Tenkai.Enzyme.API.strictAliasing!(false)
#------------------------------------------------------------------------------
xmin, xmax = 0.0, 1.0
degree = 2
γ = 5/3
final_time = 0.45 # 0.45
nx = 200
boundary_condition = (neumann, neumann)
initial_value, exact_solution = Eq.RHD_Ryu_RP_P1, Eq.exact_RHD_Ryu_RP_P1

dummy_bv(x,t) = 0.0
boundary_value = dummy_bv

solver = LWEnzymeTower()
solution_points = "gl"
correction_function = "radau" # "g2"
# To make AD work, you need Enzyme.API.strictAliasing!(false)

numerical_flux = Eq.rusanov_for_ad
bound_limit = "yes"
bflux = evaluate

cfl = 0.0
bounds = ([-Inf],[Inf]) # Not used in Euler
tvbM = 0.0
save_iter_interval = 0
save_time_interval = 0.0 #final_time/100 #0.0
animate = true # Factor on save_iter_interval or save_time_interval
compute_error_interval = 0

# blend parameters
indicator_model = "model1"
debug_blend =  false # false
cfl_safety_factor = 0.95 #0.95 #0.95
pure_fv = false #true
adjust_flux_arg = 1 #1=scale #2=replace
eos = 3  #1=ID #2=TM #3=RC #4=IP

equation = Eq.get_equation(γ, eos, adjust_flux_arg)
#------------------------------------------------------------------------------
Eq.RHD_Ryu_RP_P1(x) = Eq.RHD_Ryu_RP_P1(x,equation)
Eq.exact_RHD_Ryu_RP_P1(x,t) = Eq.exact_RHD_Ryu_RP_P1(x,t,equation)

grid_size = nx
domain = [xmin, xmax]
problem = Problem(domain, initial_value, boundary_value,
                  boundary_condition, final_time, exact_solution)

limiter_blend = setup_limiter_blend(
                              blend_type = fo_blend(equation),
                              #indicating_variables = Eq.primitive_indicator!,
                              #indicating_variables = Eq.rho_p_indicator!,
                              indicating_variables = Eq.rho_lorentz_p_indicator!, #
                              #indicating_variables = Eq.conservative_indicator!,
                              reconstruction_variables = conservative_reconstruction,
                              indicator_model = indicator_model,
                              constant_node_factor = 1.0,
                              amax = 1.0,
                              debug_blend = debug_blend,
                              pure_fv = pure_fv
                            )

limiter_none = setup_limiter_none()
#  limiter = setup_limiter_tvb(equation; tvbM = tvbM)
# limiter = setup_limiter_hierarchical(alpha = 1.0,
#                                      reconstruction = characteristic_reconstruction)
limiter = limiter_blend
scheme = Scheme(solver, degree, solution_points, correction_function,
                numerical_flux, bound_limit, limiter, bflux)
param = Tenkai.Parameters(
                    grid_size, cfl, bounds, save_iter_interval,
                    save_time_interval, compute_error_interval;
                    animate = animate,
                    cfl_safety_factor = cfl_safety_factor,
                    saveto = "none"
                    # time_scheme = "SSPRK54"
                  )
#------------------------------------------------------------------------------
# problem, scheme, param = ParseCommandLine(problem, param, scheme,
#                                           equation, ARGS) ##NN
#------------------------------------------------------------------------------
@time sol = Tenkai.solve(equation, problem, scheme, param);

#println(sol["errors"])

#return sol;
sol["plot_data"].p_ua
