
### Staring with julia

### run this line to be able to install packages
#using Pkg

### run this to install packages if necessary
#Pkg.add("DifferentialEquations")
#Pkg.add("Plots")

using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)


## Implementation of the toy example of the slides 8.4.2022
## implement it for ODE (ordinary differential equations) solver

#
function toy_SIR(dx,x,p,t)
     # dx is the time derivative x - vector valued function of time t
     # x vector valued function of individuals (nice, nasty)
     # p vector of parameters - alpha

     A, B = x
     alpha = p

     dx[1] = - alpha * A    # dA/dt
     dx[2] = alpha * A      # dB/dt
     dx

end

## Use ODE solver

## formulate inital-value problem
x0 = [1000., 0]  ## Initial condition A(0)=10000, B(0)=0  t=0 is inital starting point
tspan = (0.0 , 50.)   ## We want to solve the ODE from time t=0 to time t=200
p = 1/3  ## rate alpha alpha=1/3 days, i.e., 1/alpha =3 days

problem = ODEProblem(toy_SIR,x0,tspan,p) ## problem is a Julia Object created
# by the function ODEProblem - this object formulates the ODE proble

solution = solve(problem)  ## function solve, solves the ODE problem
### solution is a Julia Object

solution.t
solution.u


plot(solution)


### Home exercise - try to implement the Lorenz attractor
### ODE help 
