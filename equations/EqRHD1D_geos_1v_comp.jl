module EqRHD1D_geos_1v_comp

# This code has been authored by Sujoy Basak, IIT Delhi
# https://scholar.google.com/citations?user=lN1j2ggAAAAJ&hl=en

#import GZip
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
using Roots

using Tenkai
using Tenkai.Basis

import Tenkai: admissibility_tolerance

( # Methods to be extended
import Tenkai: flux, prim2con, prim2con!, con2prim, con2prim!,
             eigmatrix,
             limit_slope, zhang_shu_flux_fix,
             apply_tvb_limiter!, apply_bound_limiter!, initialize_plot,
             write_soln!, compute_time_step, post_process_soln,
            # rusanov, primitive_indicator!, rho_indicator!, rho_p_indicator, p_indicator,
             is_admissible
)

(
using Tenkai: PlotData, data_dir, get_filename, neumann, minmod,
               get_node_vars,
               set_node_vars!,
               nvariables, eachvariable,
               add_to_node_vars!, subtract_from_node_vars!,
               multiply_add_to_node_vars!
)

using Tenkai.MuladdMacro

# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin

struct RHD1D <: AbstractEquations{1,3}
   γ::Float64
   nvar::Int64
   name::String
   initial_values::Dict{String, Function}
   numfluxes::Dict{String, Function}
   eos::Int64
   adjust_flux_arg::Int64 #should it be here?
end

#-------------------------------------------------------------------------------
# PDE Information
#-------------------------------------------------------------------------------

@inbounds @inline function flux(x, U, eq::RHD1D)
   D, m, E = U
   (ρ, v, p) = con2prim(eq, U)
    # Ignore orientation since it is always "1" in 1D
    f1 = D*v
    f2 = m*v + p
    f3 = m
   return SVector(f1, f2, f3)
end

@inbounds @inline flux(U, eq::RHD1D) = flux(0.0, U, eq)

@inline @inbounds function get_enthalpy(eq::RHD1D, prim)
    @unpack eos = eq
    Theta = prim[3]/prim[1]
    if eos == 1
      h = 1+Theta*eq.γ/(eq.γ-1)
    elseif eos == 2
      h = Theta*(5/2) + sqrt(1+(9/4)*Theta^2)
    elseif eos == 3
      h = 2*(6*Theta^2+4*Theta+1)/(3*Theta+2)
    elseif eos == 4
      h = 2*Theta + sqrt(1+4*Theta^2)
    else
      @assert false "Equation of state is not supported"
    end
    return h
end

# function converting primitive variables to PDE variables
@inbounds function prim2con(eq::RHD1D, prim) # primitive, gas constant
   ρ, vx, p = prim
   v = abs(vx)
   h = get_enthalpy(eq, prim)
   U1 = ρ/sqrt(1-v^2)
   U2 = (ρ/(1-v^2))*vx*h
   U3 = (ρ/(1-v^2))*h-p
   U = SVector(U1, U2, U3)
   return U
end

@inbounds function prim2con!(eq::RHD1D, ua)
   ρ, v, p = ua
   v = abs(v)
   h = get_enthalpy(eq, ua)
   U1 = ρ/sqrt(1-v^2)
   U2 = (ρ/(1-v^2))*v*h
   U3 = (ρ/(1-v^2))*h-p
    ua[1] = U1
    ua[2] = U2
    ua[3] = U3
   return nothing
end

@inbounds function prim2con!(eq::RHD1D, ua, ua_)
   ρ, v, p = ua
   v = abs(v)
   h = get_enthalpy(eq, ua)
   U1 = ρ/sqrt(1-v^2)
   U2 = (ρ/(1-v^2))*v*h
   U3 = (ρ/(1-v^2))*h-p
    ua_[1] = U1
    ua_[2] = U2
    ua_[3] = U3
   return nothing
end

# function converting pde variables to primitive variables
@inline @inbounds function f_and_df_id(eq::RHD1D, x::Float64, u::AbstractArray) #function of |v|
   D=u[1]
   M=u[2]
   E=u[3]
   γ = eq.γ
   abM = abs(M)

   a = 2*γ*abM*E
   gm1 = γ-1
   deno = gm1^2*(abM^2+D^2)
   if deno< 1e-13
      deno +=1e-20
   end

   a3=-(a*gm1)/deno
   a2=(γ^2*E^2+2*gm1*abM^2-gm1^2*D^2)/deno
   a1=-a/deno
   a0=abM^2/deno
   f = a0 + a1*x + a2*x^2 + a3*x^3 + x^4
   df = a1 + 2*a2*x + 3*a3*x^2 + 4*x^3
   return f, df
end

@inline @inbounds function con2prim_id(eq::RHD1D, u::AbstractArray)
   @unpack γ = eq
   D = u[1]
   M = u[2]
   E = u[3]
   abM=abs(M)
   # @show M
   if abs(abM) < 1e-12  # 1e-13 creates problem with lower bound of |v|
      abv = 0.0
   else
      sss=γ^2*E^2-4*(γ-1)*abM^2
      if sss<0
         ulb = 0
      else
         ulb= (1/(2*abM*(γ-1))) * (γ*E-sqrt(sss))
      end

      deno = E
      if abs(deno) < 1e-13
         deno += 1e-20
      end
      uub= min(1,abM/deno + 1e-6 )

      # if uub < 1e-13
      #    ab_v = 0.0
      # else

      if ulb > 1e-9
         z = (1/2) * (1-D/(deno))*(ulb-uub)
      else
         z = 0
      end
      # @show ulb, uub
      abv0 = (1/2)*(ulb+uub) + z
      # @show abv0
      abv = abv0
      i = 0                     #Starting counter for NR loop
      # @show u, "second"
      f, fd = f_and_df_id(eq, abv, u)
      # dd = 1.0
      dd = -f/fd
      tol = 1e-17 #17 #13 #13 #15 #20    #Think why 1e-16 starts creating problem and decreasing it giving better result for Ryu_RP_P1
                         # < 1e-15 giving error for smooth with vx=0.99
      max_iter = 40  #Sometime needs more than 20 iterations to converge with tol 1e-20
      while abs(f) > tol && i < max_iter && abs(dd) > 1e-13 && abs(fd) > tol
         # dd = -f/fd
         abv = abv + dd    # update of Newtonraphson method
         # @show abv, f
         # abv_new = abv + dd    # update of Newtonraphson method
         # if abs(abv_new) < 1.0 # Think why is it needed
         #    abv = abv_new
         # else
         #    @show "break is running"
         #    break
         # end
         i += 1
         f, fd =  f_and_df_id(eq, abv, u)
         dd = -f/fd
      end
      if abs(f) > 1e-13 && abs(dd) > 1e-13
         println("NR method in con2prim_id did not coverge and i, f are ", i," ", f)
      end
      # ANALYTIC ############################# #Analytic does not give tol to 1e-20 always
      # a = 2*γ*abM*E
      # gm1 = γ-1
      # deno = gm1^2*(abM^2+D^2)
      # if deno< 1e-13
      #    deno +=1e-20
      # end

      # a3=-(a*gm1)/deno
      # a2=(γ^2*E^2+2*gm1*abM^2-gm1^2*D^2)/deno
      # a1=-a/deno
      # a0=abM^2/deno
      # b2=-a2
      # b1=a1*a3-4*a0
      # b0=-(a1^2+a0*a3^2-4*a0*a2)
      # r=(b1*b2-3*b0)/6-(b2^3)/27
      # q=b1/3-b2^2/9
      # # s2=cbrt(r-sqrt(abs(q^3+r^2)))
      # # s1=cbrt(r+sqrt(abs(q^3+r^2)))
      # # z1=(s1+s2)-b2/3
      # # Cd=z1/2-sqrt(abs((z1/2)^2-a0))
      # # Bd=a3/2+sqrt(abs((a3^2)/4+z1-a2))
      # # ab_v=(-Bd+sqrt(abs(Bd^2-4*Cd)))/2
      # q3r2 = q^3+r^2
      # if abs(q3r2)<1e-13
      #    q3r2 = 0.0
      # end
      # if q3r2 < -1e-6
      #    @show q3r2
      # end
      # s2=cbrt(r-sqrt(abs(q3r2)))
      # s1=cbrt(r+sqrt(abs(q3r2)))
      # z1=(s1+s2)-b2/3
      # arg1 = (z1/2)^2-a0
      # if abs(arg1) < 1e-13
      #    arg1 = 0.0
      # end
      # if arg1 < -1e-6
      #    @show arg1
      # end
      # Cd=z1/2-sqrt(abs(arg1))
      # arg2 = (a3^2)/4+z1-a2
      # if abs(arg2) < 1e-13
      #    arg2 = 0.0
      # end
      # if arg2 < -1e-6
      #    @show arg2
      # end
      # Bd=a3/2+sqrt(abs(arg2))
      # arg3 = Bd^2-4*Cd
      # if abs(arg3) < 1e-13
      #    arg3 = 0.0
      # end
      # abv=(-Bd+sqrt(abs(arg3)))/2

      # if abv <1e-15 #-10 #for BLAST
      #    abv=0.0
      # end
   end

   #To find velocity
   if abM <1e-13 #1e-15
      v = abv
   else
      v = (M/abM)*abv
   end

   #To find density
   if abv > 1.0 && abv < 1.0 + 1e-6 #Temorary fix for WU_RP1
      abv = 1.0
   end
   # ρ = D*sqrt(abs(1-abv^2)) #SEE
    ρ = D*sqrt(1-abv^2) #SEE

   #To find pressure
   mv = M*v
   p = (γ-1.0)*(E-mv-ρ)

   return SVector(ρ,v,p)
end

@inline @inbounds function con2prim_tm(eq::RHD1D, u::AbstractArray)
   M = u[2]
   abM=abs(M)
   D=u[1]
   E=u[3]
   deno=(E^2-abM^2)^2
   # if deno< 1e-13 #Is it needed?
   #    deno +=1e-20
   # end

   em = E^2 + abM^2
   md = abM^2 + D^2
   me = abM^2 * E^2

   c1 = (em*(4*em-md)-14*me)/(2*deno)
   c2 = ((4*em-md)^2-57*me)/(16*deno)
   c3 = - 9*me/(16*deno)
   J = (3*c2-c1^2)/9
   H = (9*c1*c2-27*c3-2*c1^3)/54
   cosarg = H/sqrt(-J^3)

   tt = acos(cosarg + 0.0im)
   costerm = cos(tt/3)
   if imag(costerm) > 1e-10
      @assert false "imaginary part of cos_term is large"
      @show imag(costerm)
   end
   costerm = real(costerm)
   W = 2*sqrt(-J)*costerm-c1/3
   argd = W/(W+1)
   if abs(argd) < 1e-15 #Think why becoming negative and 1e-13 is bad
      argd = 0.0
   end
   ab_v = sqrt(argd)

   #velocity
   if abM <1e-13
      v = ab_v
   else
      v = (M/abM)*ab_v
   end

   #density
   # ρ = D*sqrt(abs(1-ab_v^2))
   ρ = D*sqrt(1-ab_v^2) #SEE

   #pressure
   mv = M*v
   emv = E-mv
   if abs(emv) < 1e-13 #Is it needed?
      emv += 1e-20
   end
   p = emv/3 - ρ^2/(3*emv) #(emv^2-ρ^2)/(3*emv)

   return SVector(ρ,v,p)
end

@inline @inbounds function f_and_df_rc(x, u::AbstractArray) #function of PI
   M=abs(u[2])
   D=u[1]
   E=u[3]
   a = sqrt(x^2 - M^2)
   h = a/D
   dhdx = x/(a*D)
   b = sqrt((3*h + 8)^2 - 96)
   f = h/8 - 1/3 + b / 24
   dfdh = 1/8 + 3*(3*h + 8)/(24*b)
   F = x^2 - x * E - D^2 * h * f
   dFdx = 2 * x - E - D^2 * dhdx * (f + h * dfdh)
   return F, dFdx
end

@inline @inbounds function con2prim_rc(eq::RHD1D, u::AbstractArray)
   M=u[2]
   D=u[1]
   E=u[3]

   x = E
   i = 0
   f, fd = f_and_df_rc(x, u)
   reldx = 1.0
   max_iter = 30 #think if 10 is sufficient
   tol = 1.0e-13
   while abs(f) > tol && i < max_iter && reldx > tol
      dx = -f/fd
      reldx = abs(dx/x)
      x += dx
      i += 1
      f, fd =  f_and_df_rc(x, u)
      # if i > 20
      #    @show i
      # end
      #@printf("Iter,x,dx,dx/x,f(x) = %4d %20.18e %20.12e %20.12e %20.18e\n",i,x,dx,reldx,f)
   end
   if abs(f) > tol && reldx > tol
      @assert false "The NR method for conservative to primitive conversion did not converge"
   end

   #velocity
   v = M / x
   #pressure
   p = x - E
   #density
   rho = D * sqrt(1 - v^2)
   return SVector(rho,v,p)
end

@inline @inbounds function con2prim_ip(eq::RHD1D, u::AbstractArray)
   D = u[1]
   abM = abs(u[2])
   E = u[3]
   p = 2/3 * (sqrt(E^2 - 3/4 * (abM^2 + D^2))- 0.5* E)

   #velocity
   v = u[2]/(E + p)
   #density
   rho = D * sqrt(1 - v^2)
   return SVector(rho,v,p)
end

@inbounds @inline function con2prim(eq::RHD1D, U)
   if eq.eos == 1
      primitives = con2prim_id(eq,U)
   elseif eq.eos == 2
      primitives = con2prim_tm(eq,U)
   elseif eq.eos == 3
      primitives = con2prim_rc(eq,U)
   elseif eq.eos == 4
      primitives = con2prim_ip(eq,U)
   else
      @assert false "Equation of state is not supported"
   end
   return primitives
end

@inbounds function con2prim!(eq::RHD1D, ua, ua_)
   if eq.eos == 1
      primitives = con2prim_id(eq,ua)
   elseif eq.eos == 2
      primitives = con2prim_tm(eq,ua)
   elseif eq.eos == 3
      primitives = con2prim_rc(eq,ua)
   elseif eq.eos == 4
      primitives = con2prim_ip(eq,ua)
   else
      @assert false "Equation of state is not supported"
   end
   ua_[1] = primitives[1]
   ua_[2] = primitives[2]
   ua_[3] = primitives[3]
    return nothing
 end

 @inbounds function con2prim!(eq::RHD1D, ua)
   if eq.eos == 1
      primitives = con2prim_id(eq,ua)
   elseif eq.eos == 2
      primitives = con2prim_tm(eq,ua)
   elseif eq.eos == 3
      primitives = con2prim_rc(eq,ua)
   elseif eq.eos == 4
      primitives = con2prim_ip(eq,ua)
   else
      @assert false "Equation of state is not supported"
   end
   ua[1] = primitives[1]
   ua[2] = primitives[2]
   ua[3] = primitives[3]
   return nothing
 end

 @inline @inbounds function get_sound_speed(eq::RHD1D, ρ::Float64, v::Float64, p::Float64)
   @unpack eos = eq
   pbr = p/ρ
   if eos == 1
      c = sqrt(eq.γ*pbr*(eq.γ-1)/(eq.γ*pbr+eq.γ-1))
   elseif eos == 2
      c = sqrt((5*pbr*sqrt(pbr^2+(4/9))+3*pbr^2)/(12*pbr*sqrt(pbr^2+(4/9))+12*pbr^2+2))
   elseif eos == 3
      c = sqrt(pbr*(3*pbr+2)*(18*pbr^2+24*pbr+5)/(3*(6*pbr^2+4*pbr+1)*(9*pbr^2+12*pbr+2)))
   elseif eos == 4
      c = sqrt((2*p*sqrt(1+4*pbr^2)*ρ)/(4*p^2+4*p*sqrt(1+4*pbr^2)*ρ+ρ^2))
   else
      @assert false "Equation of state is not supported"
   end
   return c
end

@inline @inbounds function get_polytropic_index(eq::RHD1D, ρ::Float64, v::Float64, p::Float64)
   @unpack eos = eq
   pbr = p/ρ
   if eos == 1
      n = 1/(eq.γ-1)
   elseif eos == 2
      n = 3/2 + (3*pbr)/(2*sqrt(pbr^2+(4/9)))
   elseif eos == 3
      n = 3*(9*pbr^2+12*pbr+2)/(3*pbr+2)^2
   elseif eos == 4
      n = 1 + (4*p*ρ*sqrt(1+4*pbr^2))/(4*p^2+ρ^2)
   else
      @assert false "Equation of state is not supported"
   end
   return n
end

function eigmatrix(eq::RHD1D, u)
   @unpack γ = eq
   rho, v, p = con2prim(eq, u)

   h   = get_enthalpy(eq,prim)
   c   = get_sound_speed(eq, rho, v, p)
   lf  = 1/sqrt(1-v^2)
   n   = get_polytropic_index(eq, rho, v, p)

   den = h*n*c^2
   # Inverse eigenvector-matrix
   L11 = -h*(1-n*c^2)/(2*den)
   L21 = -lf*(v+n*c)/(2*den)
   L31 = lf*(1+n*c*v)/(2*den)

   L12 = h/den
   L22 = lf*v/den
   L32 = -lf/den

   L13 = -h*(1-n*c^2)/(2*den)
   L23 = -lf*(v-n*c)/(2*den)
   L33 = lf*(1-n*c*v)/(2*den)

   # L = SMatrix{nvariables(eq), nvariables(eq)}(L11, L21, L31,
   #                                             L12, L22, L32,
   #                                             L13, L23, L33)

   L = SMatrix{nvariables(eq), nvariables(eq)}(L11, L21, L31,
                                               L12, L22, L32,
                                               L13, L23, L33)'

   # Eigenvector matrix
   R11 = 1.0
   R21 = lf*h*(v-c)
   R31 = lf*h*(1-c*v)

   R12 = 1.0
   R22 = lf*h*v*(1-n*c^2)
   R32 = lf*h*(1-n*c^2)

   R13 = 1.0
   R23 = lf*h*(v+c)
   R33 = lf*h*(1+c*v)

   R = SMatrix{nvariables(eq), nvariables(eq)}(R11, R21, R31,
                                               R12, R22, R32,
                                               R13, R23, R33)
   # R = SMatrix{nvariables(eq), nvariables(eq)}(1.0, 0.0, 0.0,
   #                                             0.0, 1.0, 0.0,
   #                                             0.0, 0.0, 1.0)

   return R, L
end


#-------------------------------------------------------------------------------
# Scheme information
#-------------------------------------------------------------------------------

@inline function max_abs_eigvalue(eq::RHD1D, ρ::Float64, v::Float64, p::Float64)
    @unpack eos = eq
    c = get_sound_speed(eq, ρ, v, p)

    eig1d=(v - c)/(1 - c * v)
    eig1=abs(eig1d)
    eig2=abs(v)
    eig3d=(v + c)/(1 + c * v)
    eig3=abs(eig3d)
    eig=max(eig1,eig2,eig3)
    return eig
 end

function compute_time_step(eq::RHD1D, problem, grid, aux, op, cfl, u1, ua)
   nx = grid.size
   dx = grid.dx
   den = 0.0
   for i=1:nx
      u = get_node_vars(ua, eq, i)
      rho, v, p = con2prim(eq, u)
      smax = max_abs_eigvalue(eq,rho,v,p)
      den = max(den, smax/dx[i])
   end
   dt = cfl / den
   return dt
end

#-------------------------------------------------------------------------------
# Initial Values
#-------------------------------------------------------------------------------

function RHD_Ryu_RP_P1(x, equation)
   if x < 0.5
      rho = 10.0
      vx  = 0.0
      vy  = 0.0
      p   = 13.3
   else
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = 10.0^(-6)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Ryu_RP_P1(x, t, equation) = RHD_Ryu_RP_P1(x, equation)

function RHD_Ryu_RP_P2(x, equation)
   if x < 0.5
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = 10.0^3
   else
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = 10.0^(-2)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Ryu_RP_P2(x, t, equation) = RHD_Ryu_RP_P2(x, equation)

function RHD_Wu_RP1(x, equation)
   if x < 0.5
      rho = 1.0
      vx  = 0.0
      p   = 10.0^4
   else
      rho = 1.0
      vx  = 0.0
      p   = 10.0^(-8)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Wu_RP1(x, t, equation) = RHD_Wu_RP1(x, equation)

function RHD_Wu_shock_heating(x, equation) #DO###############
   if x < 0.5
      rho = 1.0
      vx  = 0.0
      p   = 10.0^4
   else
      rho = 1.0
      vx  = 0.0
      p   = 10.0^(-8)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Wu_shock_heating(x, t, equation) = RHD_Wu_shock_heating(x, equation)

function RHD_Wu_blast(x, equation)
   gamma = 1.4
   @unpack γ, eos = equation
   if eos == 1 && γ != gamma
      @assert false "Value of γ is not right"
   end
   if x < 0.1
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = 1000.0
   elseif x >= 0.1 && x < 0.9
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = 0.01
   else
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = 100.0
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Wu_blast(x,t,equation) = RHD_Wu_blast(x, equation)

function RHD_Migone_P1(x, equation)
   if x < 0.5
      rho = 10.0
      vx  = 0.0
      p   = 40/3
   else
      rho = 1.0
      vx  = 0.0
      p   = 2/3 * 10.0^(-6)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Migone_P1(x, t, equation) = RHD_Migone_P1(x, equation)

function RHD_Migone_P2(x, equation)
   if x < 0.5
      rho = 1.0
      vx  = 0.0
      p   = 10.0^3
   else
      rho = 1.0
      vx  = 0.0
      p   = 10.0^(-2)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Migone_P2(x, t, equation) = RHD_Migone_P2(x, equation)

function RHD_Sokolov_RP1(x, equation)
   if x < 0.5
      rho = 1.0
      vx  = 0.0
      p   = 2.0
   else
      rho = 1.0
      vx  = 0.0
      p   = 1.0
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Sokolov_RP1(x, t, equation) = RHD_Sokolov_RP1(x, equation)

function RHD_Sokolov_RP2(x, equation)
   if x < 0.5
      rho = 10.0
      vx  = 0.0
      p   = 100.0
   else
      rho = 1.0
      vx  = 0.0
      p   = 1.0
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Sokolov_RP2(x, t, equation) = RHD_Sokolov_RP2(x, equation)

function RHD_Sokolov_RP3(x, equation)
   if x < 0.4
      rho = 10.0
      vx  = 0.0
      p   = 13.3
   else
      rho = 1.0
      vx  = 0.0
      p   = 10.0^(-6)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_Sokolov_RP3(x, t, equation) = RHD_Sokolov_RP3(x, equation)

###### PROBLEMS FROM OUR FIRST PAPER ##################################

function RHD_smooth(x, equation)
    gamma = 5/3
    @unpack γ, eos = equation
    if γ != gamma
      @assert false "Value of γ is not right"
    end

    rho = 1.0 + 0.99999*sin(2*pi*x)
    vx = 0.99
    vy = 0.0
    p = 0.01

    prim=[rho, vx, p]
    U = prim2con(equation, prim)
   return U
end

exact_RHD_smooth(x,t, equation) = RHD_smooth(x-0.99*t, equation)

function RHD_smooth2(x, equation)
   @unpack γ = equation
   rho = 2.0 + sin(2*pi*x)
   vx = 0.5
   vy = 0.0
   p = 1.0
   prim=[rho, vx, p]
   U = prim2con(equation, prim)
  return U
end

exact_RHD_smooth2(x,t, equation) =  RHD_smooth2(x-0.5*t, equation)

function RHD_isentropic(x, equation)
   @unpack γ = equation
   rho_ref = 1.0
   vx_ref = 0.0
   p_ref = 100.0
   α = 1.0
   L = 0.3
   f = (abs(x) < L)*((x/L)^2 - 1)^4
   rho = rho_ref * (1.0 + α * f)

   K = p_ref/(rho_ref^γ)
   p = K*rho^γ

   h = 1+(γ*p)/(rho*(γ-1))
   c = sqrt(γ*p/(rho*h))
   A = (sqrt(γ-1)+c)/(sqrt(γ-1)-c)
   h_ref = 1 + (γ*p_ref)/(rho_ref*(γ-1))
   c_ref = sqrt(γ*p_ref/(rho_ref*h_ref))
   A_ref = (sqrt(γ-1)+c_ref)/(sqrt(γ-1)-c_ref)
   S = (A/A_ref)^(2/sqrt(γ-1))
   vx = (abs(x) < L)*(S - 1)/(S + 1)# + vx_ref ##SEE

   prim=[rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

function exact_RHD_isentropic(x, t, equation)
   L = 0.3
   @unpack γ = equation
   alph = 1.0
   rho_ref = 1.0
   p_ref = 100.0
   v_ref = 0.0
   h_ref = 1 + γ/(γ-1) * p_ref/rho_ref
   c_ref = sqrt(γ*p_ref/(rho_ref*h_ref))

   f(y) = (abs(y)<L)*((y/L)^2 -1)^4

   rho0(y)= rho_ref * (1+ alph* f(y))
   K = p_ref/(rho_ref^γ)
   p0(y) = K* rho0(y)^γ

   h0(y) = 1 + γ/(γ-1) * p0(y)/rho0(y)

   c0(y) = sqrt(γ*p0(y)/(rho0(y)*h0(y)))

   k1 = sqrt(γ-1)

   A(y) = (k1 + c0(y))/(k1 - c0(y))
   A_ref = (k1 + c_ref)/(k1 - c_ref)

   v0(y) = (abs(y)<L) * ((A_ref/A(y))^(-2/k1) - 1)/((A_ref/A(y))^(-2/k1) + 1) +v_ref ##SEE

   foo(y) = y - x + t*(v0(y)+ c0(y))/(1 + v0(y)*c0(y))
   xi = find_zero(foo,  x, atol=1e-13, rtol=1e-13)

   ##############################
   if abs(foo(xi)) > 1e-10
      @show foo(xi)
   end
   ##############################

   Jminus = 0.5*log((1+v0(x)) ./(1-v0(x))) -1/k1.*log(A(x))
   Jplus = 0.5*log((1+v0(xi)) ./(1-v0(xi)))+1/k1.*log(A(xi))

   v = (exp(Jplus+Jminus)-1) ./(exp(Jplus+Jminus)+1)
   cs = k1*(exp(0.5*k1*(Jplus-Jminus))-1) ./(exp(0.5*k1*(Jplus-Jminus))+1)
   rho = (K*(γ ./(cs .^2) - γ/(γ-1))) .^(1/(1-γ))
   p = K*rho .^γ

   prim=[rho, v, p] ##SEE
   U = prim2con(equation, prim)

   return U
end

function RHD_riemann1(x, equation)
   # gamma = 5/3
   @unpack γ = equation
   # if γ != gamma
   #    @assert false "Value of γ is not right"
   # end
   if x < 0.5
      rho = 1.0
      vx  = -0.6
      vy  = 0.0
      p   = 10.0
   else
      rho = 10.0
      vx  = 0.5
      vy  = 0.0
      p   = 20.0
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_riemann1(x, t, equation)= RHD_riemann1(x, equation)

function RHD_riemann2(x, equation)
   gamma = 5/3
   @unpack γ, eos = equation
   if γ != gamma
      @assert false "Value of γ is not right"
   end
   if x < 0.5
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = (10.0)^3
   else
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = (10.0)^(-2)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_riemann2(x, t, equation)= RHD_riemann2(x, equation)


function RHD_riemann3(x, equation)
   gamma = 5/3
   @unpack γ, eos = equation
   # if γ != gamma
   #    @assert false "Value of γ is not right"
   # end
   if x < 0.5
      rho = 10.0
      vx  = 0.0
      vy  = 0.0
      p   = 40/3
   else
      rho = 1.0
      vx  = 0.0
      vy  = 0.0
      p   = (10.0)^(-6)   #2/3 * (10.0)^(-5)
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_riemann3(x, t, equation)= RHD_riemann3(x, equation)

function RHD_density_pert(x, equation)
   # gamma = 5/3
   @unpack γ, eos = equation
   # if γ != gamma
   #    @assert false "Value of γ is not right"
   # end
   if x < 0.5
      rho = 5.0
      vx  = 0.0
      vy  = 0.0
      p   = 50.0
   else
      rho = 2.0 + 0.3 * sin(50 * x)
      vx  = 0.0
      vy  = 0.0
      p   = 5.0
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_density_pert(x,t, equation)= RHD_density_pert(x, equation)

function RHD_riemann4(x, equation)
   gamma = 4/3
   @unpack γ, eos = equation
   if γ != gamma
      @assert false "Value of γ is not right"
   end
   if x < 0.5
      rho = 1.0
      vx  = 0.9
      p   = 1.0
   else
      rho = 1.0
      vx  = 0.0
      p   = 10.0
   end
   prim = [rho, vx, p]
   U = prim2con(equation, prim)
   return U
end

exact_RHD_riemann4(x, t, equation)= RHD_riemann4(x, equation)
###############################################################



initial_values = Dict{String, Function}()
initial_values["RHDsmooth"] = RHD_smooth
initial_values["RHDsmooth2"] = RHD_smooth2
initial_values["RHD_isentropic"] = RHD_isentropic
initial_values["RHD_RP1"] = RHD_riemann1
initial_values["RHD_RP2"] = RHD_riemann2
initial_values["RHD_RP3"] = RHD_riemann3
initial_values["RHD_density_pert"] = RHD_density_pert
initial_values["RHD_RP4"] = RHD_riemann4
# initial_values["RHD_blast"] = RHD_blast
initial_values["RHD_Ryu_RP_P1"] = RHD_Ryu_RP_P1
initial_values["RHD_Ryu_RP_P2"] = RHD_Ryu_RP_P2
initial_values["RHD_Wu_RP1"] = RHD_Wu_RP1
initial_values["RHD_Wu_shock_heating"] = RHD_Wu_shock_heating
initial_values["RHD_Wu_blast"] = RHD_Wu_blast
initial_values["RHD_Migone_P1"] = RHD_Migone_P1
initial_values["RHD_Migone_P2"] = RHD_Migone_P2
initial_values["RHD_Sokolov_RP1"] = RHD_Sokolov_RP1
initial_values["RHD_Sokolov_RP2"] = RHD_Sokolov_RP2
initial_values["RHD_Sokolov_RP3"] = RHD_Sokolov_RP3
#-------------------------------------------------------------------------------
# Numerical Fluxes
#-------------------------------------------------------------------------------



@inbounds @inline function rusanov(x, ual, uar, Fl, Fr, Ul, Ur, eq::RHD1D, dir)
    rho_ll,v_ll,p_ll=con2prim(eq,ual)
    rho_rr,v_rr,p_rr=con2prim(eq,uar)

   λ = max(max_abs_eigvalue(eq,rho_ll,v_ll,p_ll),max_abs_eigvalue(eq,rho_rr,v_rr,p_rr)) # local wave speed

   return 0.5 * (Fl + Fr - λ * (Ur - Ul))
end

function rusanov_for_ad(x, ual, uar, Fl, Fr, Ul, Ur, eq::RHD1D, dir)
    rusanov(x, ual, uar, Fl, Fr, Ul, Ur, eq, dir)
end

#-------------------------------------------------------------------------------
# Limiters
#-------------------------------------------------------------------------------
function apply_bound_limiter!(eq::RHD1D, grid, scheme, param, op, ua,
                                   u1, aux)
   if scheme.bound_limit == "no"
      return nothing
   end

   @timeit aux.timer "Positivity limiter" begin
   @unpack Vl, Vr = op
   nx   = grid.size
   nd   = op.degree + 1

   # Looping over tuple of functions like this is only type stable if
   # there are only two. For tuples with than two functions, see
   # https://github.com/trixi-framework/Trixi.jl/blob/0fd86e4bd856d894de6a7514edcb9758bf6f8e1e/src/callbacks_stage/positivity_zhang_shu.jl#L39

   #For D
   eps=1e-13
   for i=1:nx
      ua_=get_node_vars(ua,eq,i)
      eps=min(eps,ua_[1])
   end
   if eps < 0.0
      println("Fatal: Negative states in cell averages")
      println("D minimum cell average = $eps") ###############
      throw(DomainError(eps, "Positivity limiter failed"))
   end
   #FOR q
   for i=1:nx
      ua_=get_node_vars(ua,eq,i)
      m  = abs(ua_[2])
      q=ua_[3]-sqrt(ua_[1]^2+m^2)
       if q < 0.0
          @show i, u1[:,:,i]
          @show i, ua[:,i]
       end
      eps=min(eps,q)
   end
   if eps < 0.0
      println("Fatal: Negative states in cell averages")
      println("q minimum cell average = $eps") ###############
      throw(DomainError(eps, "Positivity limiter failed"))
   end

   #For D
   for element in 1:nx
      var_ll=var_rr=0.0
      var_min=1e20
      for i in Base.OneTo(nd)
         u_node= get_node_vars(u1, eq, i, element)
         var=u_node[1]
         var_ll += var * Vl[i]
         var_rr += var * Vr[i]
         var_min = min(var_min, var)
      end
      var_min= min(var_min, var_ll, var_rr)
      ua_ = get_node_vars(ua, eq, element)
      var_avg=ua_[1]
      dd1=abs(var_min - var_avg)
      if dd1< 1e-13
         dd1 += 1e-20
      end
      ratio = abs(eps - var_avg)/dd1
      theta = min(ratio, 1.0)
      if theta < 1.0
         for i=1:nd
            u_node = get_node_vars(u1, eq, i, element)
            multiply_add_set_node_vars!(u1,
                                          theta, u_node,
                                          1-theta, ua_,
                                          eq, i, element)
         end
      end
   end
   #For q
   for element in 1:nx
      var_ll=var_rr=0.0
      udl = zeros(3)######
      udr = zeros(3)######
      var_min=1e20
      for i in Base.OneTo(nd)
         u_node= get_node_vars(u1, eq, i, element)
         mm =  abs(u_node[2])
         var=u_node[3]-sqrt(abs(u_node[1]^2+mm^2))
         var_ll += var * Vl[i] ##
         var_rr += var * Vr[i] ##
         udl += u_node .* Vl[i] #######
         udr += u_node .* Vr[i] #######
         var_min = min(var_min, var)
      end
      var_ll_s = udl[3]-sqrt(udl[1]^2+udl[2]^2) ######
      var_rr_s = udr[3]-sqrt(udr[1]^2+udr[2]^2) ######
      var_min= min(var_min, var_ll, var_rr, var_ll_s, var_rr_s) ###### TEMPORARY SEE
      #println("var_ll= $var_ll, var_ll_s= $var_ll_s")
      ua_ = get_node_vars(ua, eq, element)
      mmd = abs(ua_[2])
      var_avg=ua_[3]-sqrt(ua_[1]^2 + mmd^2)
      #var_avg = ua_[1]
      dd2=abs(var_min - var_avg)
      if dd2< 1e-13
         dd2 += 1e-20
      end
      #ratio = abs(eps - var_avg)/(abs(var_min - var_avg) + 1e-13)
      ratio = abs(eps - var_avg)/dd2
      theta = min(ratio, 1.0)
      if theta < 1.0
         for i=1:nd
            u_node = get_node_vars(u1, eq, i, element)
            multiply_add_set_node_vars!(u1,
                                          theta, u_node,
                                          1-theta, ua_,
                                          eq, i, element)
         end
      end
   end
   end #timer
end

function apply_tvb_limiter!(eq::RHD1D, problem, scheme, grid, param, op, ua,
                                 u1, aux)
   @timeit aux.timer "TVB limiter" begin
   nx = grid.size
   @unpack xg, wg, Vl, Vr = op
   @unpack limiter = scheme
   @unpack tvbM, cache = limiter
   left_bc, right_bc = problem.boundary_condition
   nd = length(wg)
   nvar = nvariables(eq)
   # face values
   (uimh, uiph ,Δul,Δur,Δual,Δuar,char_Δul,char_Δur,char_Δual,char_Δuar,
    dulm,durm, du ) = cache


   # Loop over cells
   for cell in 1:nx
      ual, ua_, uar = (get_node_vars(ua, eq, cell-1),
                       get_node_vars(ua, eq, cell),
                       get_node_vars(ua, eq, cell+1))
      R, L = eigmatrix(eq, ua_)
      fill!(uimh, zero(eltype(uimh))); fill!(uiph, zero(eltype(uiph)))
      Mdx2 = tvbM * grid.dx[cell]^2
      if left_bc == neumann && right_bc == neumann && (cell == 1 || cell == nx)
         Mdx2 = 0.0 # Force TVD on boundary for Shu-Osher ###WHY
      end
      # end # timer
      for ii=1:nd
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

      set_node_vars!(Δul ,  ua_ , eq, 1)
      set_node_vars!(Δur , uiph_, eq, 1)
      set_node_vars!(Δual,  ua_ , eq, 1)
      set_node_vars!(Δuar,  uar , eq, 1)

      subtract_from_node_vars!(Δul,  uimh_, eq) ##HOW IS IT RUNNING WITHOUT INDICES
      subtract_from_node_vars!(Δur,  ua_   , eq)
      subtract_from_node_vars!(Δual, ual   , eq)
      subtract_from_node_vars!(Δuar, ua_   , eq)

      Δul_  = get_node_vars(Δul, eq, 1)
      Δur_  = get_node_vars(Δur, eq, 1)
      Δual_ = get_node_vars(Δual, eq, 1)
      Δuar_ = get_node_vars(Δuar, eq, 1)
      mul!(char_Δul, L, Δul_)   # char_Δul = L*Δul
      mul!(char_Δur, L, Δur_)   # char_Δur = L*Δur
      mul!(char_Δual, L, Δual_) # char_Δual = L*Δual
      mul!(char_Δuar, L, Δuar_) # char_Δuar = L*Δuar

      char_Δul_ = get_node_vars(char_Δul, eq, 1)
      char_Δur_ = get_node_vars(char_Δur, eq, 1)
      char_Δual_ = get_node_vars(char_Δual, eq, 1)
      char_Δuar_ = get_node_vars(char_Δuar, eq, 1)
      for n in eachvariable(eq)
         dulm[n] = minmod(char_Δul_[n], char_Δual_[n], char_Δuar_[n], Mdx2) ##SEE HOW
         durm[n] = minmod(char_Δur_[n], char_Δual_[n], char_Δuar_[n], Mdx2)
      end

      # limit if jumps are detected
      dulm_ = get_node_vars(dulm, eq, 1)
      durm_ = get_node_vars(durm, eq, 1)
      jump_l = jump_r = 0.0
      for n=1:nvar
         jump_l += abs(char_Δul_[n]-dulm_[n])
         jump_r += abs(char_Δur_[n]-durm_[n])
      end
      jump_l /= nvar
      jump_r /= nvar

      if jump_l > 1e-06 || jump_r > 1e-06
         add_to_node_vars!(durm, dulm_, eq, 1) # durm = durm + dulm
         # We want durm = 0.5 * (dul + dur), we adjust 0.5 later
         mul!(du, R, durm)            # du = R * (dulm+durm)
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
# Blending Limiter
#-------------------------------------------------------------------------------
@inbounds @inline function primitive_indicator!(un, eq::RHD1D)
   @views for ix=1:size(un,2) # loop over dofs and faces
      @views prim = con2prim(eq,un[:,ix])
      un[1,ix] = prim[1]
      un[2,ix] = prim[2]
      un[3,ix] = prim[3]
   end
   n_ind_var = nvariables(eq)
   return n_ind_var
end

@inbounds @inline function rho_indicator!(ue, eq::RHD1D)
   # Use only density as indicating variables
   for ix=1:size(un,2) # loop over dofs and faces
      @views prim = con2prim(eq,un[:,ix])
      rho = prim[1]
      un[1,ix] = rho
   end
   n_ind_var = 1
   return n_ind_var
end

@inbounds @inline function rho_p_indicator!(un, eq::RHD1D)
   for ix=1:size(un,2) # loop over dofs and faces
      @views prim = con2prim(eq,un[:,ix])
      ρ=prim[1]
      p=prim[3]
    un[1,ix]=ρ*p
   end
   n_ind_var = 1
   return n_ind_var
end

@inbounds @inline function p_indicator!(un, eq::RHD1D)
   for ix=1:size(un,2) # loop over dofs and faces
      @views prim = con2prim(eq,un[:,ix])
      p = prim[3]
      un[1,ix] = p
   end
   n_ind_var = 1
   return n_ind_var
end

@inbounds @inline function rho_lorentz_indicator!(un, eq::RHD1D)
   for ix=1:size(un,2) #ix=2:(size(un,2)-1) # # loop over dofs and faces
      @views prim = con2prim(eq,un[:,ix])
      ρ=prim[1]
      v = prim[2]
    lf = 1/sqrt(1-v^2)
    un[1,ix]=ρ*lf
   end
   n_ind_var = 1
   return n_ind_var
end

@inbounds @inline function rho_lorentz_p_indicator!(un, eq::RHD1D)

   for ix=1:size(un,2) #for ix=2:(size(un,2)-1) # # loop over dofs and faces
      @views prim = con2prim(eq,un[:,ix])
      ρ  = prim[1]
      v = prim[2]
      lf = 1/sqrt(1-v^2)
      p  = prim[3]
      un[1,ix]=ρ*lf*p
   end
   n_ind_var = 1
   return n_ind_var
end

function is_admissible(eq::RHD1D, u::AbstractVector)
    md = abs(u[2])
    q=u[3]-sqrt(abs(u[1]^2+md^2))
   if u[1] > 0.0 && q > 0.0
      return true
   else
      return false
   end
end

function show_admissibility(u::AbstractVector) # Used for debugging
   md = abs(u[2])
   q=u[3]-sqrt(abs(u[1]^2+md^2))
    @show u[1], q
end

admissibility_tolerance(eq::RHD1D) = 1e-13

function limit_slope(eq::RHD1D, s, ufl, u_s_l, ufr, u_s_r, ue, xl, xr)
    # The MUSCL-Hancock scheme is guaranteed to be admissibility preserving if
   # slope is chosen so that
   # u_star_l = ue + 2.0*slope*xl, u_star_r = ue+2.0*slope*xr are admissible
   # ue is already admissible and we know we can find sequences of thetas
   # to make theta*u_star_l+(1-theta)*ue is admissible.
   # This is equivalent to replacing u_star_l by
   # u_star_l = ue + 2.0*theta*s*xl.
   # Thus, we simply have to update the slope by multiplying by theta.
   #@unpack γ = eq
   eps = admissibility_tolerance(eq)

   #FOR D
   var_star_tuple= (u_s_l[1],u_s_r[1])
   var_low= ue[1]
   theta = 1.0
      for var_star in var_star_tuple
         if var_star < eps
            # TOTHINK - Replace eps here by 0.1*var_low
            ddd1=abs(var_star - var_low)
            if ddd1 < 1e-13 ##SEE
               ddd1 += 1e-20
            end
            #ratio = abs(0.1*var_low - var_low) / (abs(var_star - var_low) + 1e-13 )
            ratio = abs(0.1*var_low - var_low) / ddd1
            theta = min(ratio, theta)
         end
      end
      s *= theta
      u_s_l = ue + 2.0*xl*s
      u_s_r = ue + 2.0*xr*s

   #FOR q
   ml = abs(u_s_l[2])
   mr = abs(u_s_r[2])
   duml=u_s_l[3]-sqrt(u_s_l[1]^2+ml^2)
   dumr=u_s_r[3]-sqrt(u_s_r[1]^2+mr^2)
   var_star_tuple= (duml,dumr)
   mlow = abs(ue[2])
   var_low= ue[3]-sqrt(ue[1]^2+mlow^2)
   theta = 1.0
      for var_star in var_star_tuple
         if var_star < eps
            # TOTHINK - Replace eps here by 0.1*var_low
            ddd2=abs(var_star - var_low)
            if ddd2 < 1e-13
               ddd2 += 1e-20
            end
            ratio = abs(0.1*var_low - var_low) / ddd2
            theta = min(ratio, theta)
         end
      end
      s *= theta
      u_s_l = ue + 2.0*xl*s
      u_s_r = ue + 2.0*xr*s

   ufl = ue + xl*s
   ufr = ue + xr*s

   return ufl, ufr
end

function zhang_shu_flux_fix(eq::RHD1D,
                                 uprev,    # Solution at previous time level
                                 ulow,     # low order update
                                 Fn,       # Blended flux candidate
                                 fn_inner, # Inner part of flux
                                 fn,       # low order flux
                                 c         # c is such that unew = u - c(fr-fl)
                                 )
uhigh = uprev - c * (Fn-fn_inner) # First candidate for high order update
D_low=ulow[1]
# if D_low < 0.0
#    @show D_low
# end
D_high=uhigh[1]
eps = 0.1*D_low
dn=abs(D_high - D_low)
if dn < 1e-13
   dn += 1e-20
end
ratio = abs(eps - D_low)/dn
theta = min(ratio,1.0)
if theta < 1.0
   Fn= theta*Fn + (1.0-theta)*fn # Second candidate for flux
end
uhigh = uprev - c * (Fn - fn_inner) # Second candidate for uhigh

D_high_test = uhigh[1]
if D_high < 0.0
   @show D_low, D_high, D_high_test
end

mll = abs(ulow[2])
mhh = abs(uhigh[2])
q_low=ulow[3]- sqrt(ulow[1]^2+mll^2)

q_high=uhigh[3]- sqrt(uhigh[1]^2+mhh^2)

eps=0.1*q_low
dn2=abs(q_high-q_low)

if dn2 < 1e-13 #1e-10 ##SEE
   dn2 += 1e-20
end
ratio= abs(eps-q_low)/dn2
theta = min(ratio, 1.0)
if theta < 1.0
   Fn = theta*Fn + (1.0-theta)*fn #Final flux
end

return Fn
end


function initialize_plot(eq::RHD1D, op, grid, problem, scheme, timer, u1,
   ua)
@timeit timer "Write solution" begin
@timeit timer "Initialize write solution" begin
# Clear and re-create output directory
rm("output", force=true, recursive=true)
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
grid = false, showaxis = false, bottom_margin = 0Plots.px);
# Initialize subplots for density, velocity and pressure

p_u1 = [plot() for _=1:nvar];
p_ua = [plot() for _=1:nvar];
labels = ["Density", "Velocity", "Pressure"]




y = zeros(nx) # put dummy to fix plotly bug with OffsetArrays ##SEE

for n=1:nvar
@views plot!(p_ua[n], xc, y, label = "Approximate",
linestyle = :dot, seriestype = :scatter,
color = :blue, markerstrokestyle = :dot,
markershape = :circle, markersize = 2, markerstrokealpha = 0);
xlabel!(p_ua[n], "x"); ylabel!(p_ua[n], labels[n])
end

l_super = @layout[ a{0.01h}; b c d] # Selecting layout for p_title being title
p_ua = plot(p_title, p_ua[1], p_ua[2], p_ua[3], layout = l_super,
size = (1500,500)); # Make subplots




# Set up p_u1 to contain polynomial approximation as a different curve
# for each cell
x   = LinRange(xf[1], xf[2], nu)
up1 = zeros(nvar, nd)
u   = zeros(nu)

# for ii=1:nd
#    @views con2prim!(eq, u1[:,ii,1],up1[:,ii]) # store prim form in up1
# end

@views up1[:,:] .= u1[:,:,1]

#xgn=xf[1] .+ xg*(xf[2]-xf[1])

for n=1:nvar
   u = @views Vu * up1[n,:]
   #u = @views up1[n,:]
   plot!(p_u1[n], x, u, color = :red, legend=false, label="for nx=$nx",lw=2);
   xlabel!(p_u1[n], "x"); ylabel!(p_u1[n], labels[n]);
end
for i=2:nx
   # for ii=1:nd
   #    @views con2prim!(eq, u1[:,ii,i], up1[:,ii]) # store prim form in up1
   # end
   @views up1[:,:] .= u1[:,:,i]
   x = LinRange(xf[i], xf[i+1], nu)

   #xgn=xf[i] .+ xg*(xf[i+1]-xf[i])

   for n=1:nvar
      u = @views Vu * up1[n, :]
      #u = @views up1[n,:]
      plot!(p_u1[n], x, u, color = :red, label="for nx=$nx", legend=false,lw=2) ##dekh nothing was in label
   end
end

l = @layout[ a{0.01h}; b c d] # Selecting layout for p_title being title
p_u1 = plot!(p_title, p_u1[1], p_u1[2], p_u1[3], layout = l,
size = (1700,500)); # Make subplots ##!


anim_ua, anim_u1 = Animation(), Animation(); # Initialize animation objects
plot_data = PlotData(p_ua, anim_ua, p_u1, anim_u1);
return plot_data
end # timer
end # timer
end



function write_soln!(base_name, fcount, iter, time, dt, eq::RHD1D, grid,
    problem, param, op, ua, u1, aux, ndigits=3)
   @timeit aux.timer "Write solution" begin
   @unpack plot_data = aux
   avg_filename = get_filename("output/avg", ndigits, fcount)
   @unpack p_ua, p_u1, anim_ua, anim_u1 = plot_data
   @unpack final_time = problem
   xc = grid.xc
   nx = grid.size
   @unpack xg = op
   nd = op.degree + 1
   nu = nd
   #xu = LinRange(0.0, 1.0, nu)
   #Vu = Vandermonde_lag(xg, xu)
   xu = xg ##
   nvar = eq.nvar
   @unpack save_time_interval, save_iter_interval, animate = param
   avg_file = open("$avg_filename.txt", "w")
   up_ = zeros(nvar)
   ylims = [[Inf,-Inf] for _=1:nvar] # set ylims for plots of all variables
   for i=1:nx
      @views con2prim!(eq, ua[:,i],up_) # store primitve form in up_
      @printf(avg_file, "%e %e %e %e\n", xc[i], up_[1], up_[2], up_[3])
      # TOTHINK - Check efficiency of printf
      for n=1:eq.nvar
         p_ua[n+1][1][:y][i] = @views up_[n]    # Update y-series ## SINCE FIRST ONE IS TITLE PLOT
         ylims[n][1] = min(ylims[n][1], up_[n]) # Compute ymin
         ylims[n][2] = max(ylims[n][2], up_[n]) # Compute ymax
      end
   end
   close(avg_file)
   for n=1:nvar # set ymin, ymax for ua, u1 plots
   ylims!(p_ua[n+1],(ylims[n][1]-0.1,ylims[n][2]+0.1))
   ylims!(p_u1[n+1],(ylims[n][1]-0.1,ylims[n][2]+0.1))
   end
   #println("ylimu=$(ylims[1][2]-0.1)")
   #println("yliml=$(ylims[1][1]-0.1)")
   t = round(time; digits=3)
   title!(p_ua[1], "Cell averages plot, $nx cells, t = $t")
   sol_filename = get_filename("output/sol", ndigits, fcount)
   sol_file = open(sol_filename*".txt", "w")
   up1 = zeros(nvar,nd)

   u = zeros(nvar,nu)
   x = zeros(nu)

   for i=1:nx
      for ii=1:nd
         @views con2prim!(eq, u1[:,ii,i], up1[:,ii]) # store prim form in up1
      end
      @. x = grid.xf[i] + grid.dx[i]*xu
      #@views mul!(u, up1, Vu')
      u = up1 ##
      for n=1:nvar
         p_u1[n+1][i][:y] = u[n,:]
      end
      for ii=1:nu
         @printf(sol_file, "%e %e %e %e\n", x[ii], u[1,ii], u[2,ii], u[3,ii])
      end
   end

   close(sol_file)
   title!(p_u1[1], "Numerical Solution, $nx cells, t = $t")
   println("Wrote $sol_filename.txt, $avg_filename.txt")
   if problem.final_time - time < 1e-10
      cp("$avg_filename.txt","./output/avg.txt", force=true)
      cp("$sol_filename.txt","./output/sol.txt", force=true)
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
         k1, k2 = ceil(time/animate_time_interval), floor(time/animate_time_interval)
         if (abs(time-k1*animate_time_interval) < 1e-10 ||
            abs(time-k2*animate_time_interval) < 1e-10)
            frame(anim_ua, p_ua)
            frame(anim_u1, p_u1)
         end
      end
   end
   fcount += 1
   return fcount
   end # timer
end

function post_process_soln(eq::RHD1D, aux, problem, param, scheme)
   @unpack timer, error_file = aux
   @timeit timer "Write solution" begin
   println("Post processing solution")
   nvar = eq.nvar
   @unpack plot_data = aux
   @unpack p_ua, p_u1, anim_ua, anim_u1 = plot_data
   @unpack animate, saveto = param
   initial_values = eq.initial_values
   if problem.initial_value in values(initial_values) # Using ready made tests
      initial_value_string, = [a for (a,b) in initial_values if
                              b==problem.initial_value]
   end

   savefig(p_ua, "output/avg.png")
   if scheme.degree > 0
      savefig(p_u1, "output/sol.png")
      savefig(p_u1, "output/sol.html")
   end
   savefig(p_ua, "output/avg.html")
   if animate == true
      gif(anim_ua, "output/avg.mp4", fps = 5)
      gif(anim_u1, "output/sol.mp4", fps = 5)
   end
   println("Wrote avg, sol in gif,html,png format to output directory.")
   plot(p_ua); plot(p_u1);

   close(error_file)
   if saveto != "none"
      if saveto[end] == "/"
         saveto = saveto[1:end-1]
      end
      mkpath(saveto)
      for file in readdir("./output")
         cp("./output/$file", "$saveto/$file", force=true)
      end
      cp("./error.txt", "$saveto/error.txt", force=true)
      println("Saved output files to $saveto")
   end

   end # timer

   # Print timer data on screen
   print_timer(aux.timer, sortby = :firstexec); print("\n")
   show(aux.timer); print("\n")
   println("Time outside write_soln = "
            * "$(( TimerOutputs.tottime(timer)
                  - TimerOutputs.time(timer["Write solution"]) ) * 1e-9)s")
   println("─────────────────────────────────────────────────────────────────────────────────────────")
   timer_file = open("./output/timer.json", "w")
   JSON3.write(timer_file, TimerOutputs.todict(timer))
   close(timer_file)
   return nothing
end

function get_equation(γ, eos, adjust_flux_arg)
   name = "1d RHD Equations"
   numfluxes = Dict("rusanov"  => rusanov)
   nvar = 3
   return RHD1D(γ, nvar,
                  name, initial_values, numfluxes, eos, adjust_flux_arg)
end


end # @muladd

end