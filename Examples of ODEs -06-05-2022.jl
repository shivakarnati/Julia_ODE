### Staring with julia

### run this line to be able to install packages
#using Pkg

### run this to install packages if necessary
#Pkg.add("DifferentialEquations")
#Pkg.add("Plots")

using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)

#### Julia has a convenient syntax to make assignments componentwise

p .= [4, 5,7 ]   ##  3 element vector


a, b, c =  p   ## the values a, b and c are assigned within one line of code

### This would also work:

a = p[1]
b = p[2]
c = p[3]

###  line 19 is doing the same as lines 23-25


## Lorenz attractor 22.4.2022

## First function that defines the system:
## 1st version
function Lorenz_attractor1(dxyz,xyz,p,t)
     # dxyz is the time derivative of the vector (X,Y,Z) - vector valued function of time t
     # xyz is the vector (X,Y,Z)
     # p vector of parameters: (a,b,c)

     a, b, c = p

     dxyz[1] = a * (xyz[2] - xyz[1])                # dX/dt
     dxyz[2] = xyz[1] * (b - xyz[3]) - xyz[2]       # dY/dt
     dxyz[3] = xyz[1] * xyz[2] - c * xyz[3]         # dZ/dt
     dxyz

end

## 2nd version
function Lorenz_attractor2(dxyz,xyz,p,t)
     # dxyz is the time derivative of the vector (X,Y,Z) - vector valued function of time t
     # xyz is the vector (X,Y,Z)
     # p vector of parameters: (a,b,c)

     a, b, c = p
     X, Y, Z = xyz

     dxyz[1] = a * (Y - X)          # dX/dt
     dxyz[2] = X * (b - Z) - Y      # dY/dt
     dxyz[3] = X * Y - c* Z         # dZ/dt
     dxyz

end


## 3rd version
function Lorenz_attractor3(dxyz,xyz,p,t)
     # dxyz is the time derivative of the vector (X,Y,Z) - vector valued function of time t
     # xyz is the vector (X,Y,Z)
     # p vector of parameters: (a,b,c)

     a, b, c = p
     X, Y, Z = xyz

     dxyz .= [ a * (Y - X) ,    X * (b - Z) - Y   ,    X * Y - c* Z ]    # (dX/dt, dY/dt, dZ/dt)
     ## fused assignment - assign each element separtely  -  recomendable in Julia
     ## it is not vector-valued and saves memory
     dxyz

end

## Use ODE solver

## formulate inital-value problem
xyz0 = [0.1, 0., 0.]   ## Initial condition X(0)=0, Y(0)=0, Z(0)=0
tspan = (0.0 , 150.)   ## We want to solve the ODE from time t=0 to time t=50
p .= [10, 28, 8/3]

prob1 = ODEProblem(Lorenz_attractor1,xyz0,tspan,p) ## problem is a Julia Object created
# by the function ODEProblem - this object formulates the ODE proble

solution = solve(prob1)  ## function solve, solves the ODE problem
### solution is a Julia Object

solution.t
solution.u


plot(solution)    ##this plots the first component of solution solution.t against solution.u (2nd component of solution)
                        # the second component is vector of 3 dimensional vectors - will be plotted side by side
                        # i.e., (t, X(t)), (t,Y(t)), (t,Z(t))


plot(solution,vars=(1,2,3)) ## variable 0 is time t, variable 1, 2 & 3 are X, Y Z, here we plot
                                                # only (X,Y,Z) for all time points, no eplicit time axis
plot(solution,vars=(0,1))  ## here we plot t and X

plot(solution,vars=(0,2))  ## here we plot t and y

plot(solution,vars=(2,3))  ## here we plot Y and Z  - no explicit time axis

plot(solution,vars=(0,2,3))  ## no plot (t,Y(t),Z(t))  - explicit time axis
### Home exercise - try to implement the Lorenz attractor
### ODE help


###

xyz0 = [0.1, 0., 0.]   ## Initial condition X(0)=0, Y(0)=0, Z(0)=0
tspan = (0.0 , 150.)   ## We want to solve the ODE from time t=0 to time t=50
p = [10, 28, 8/3]
prob2 = ODEProblem(Lorenz_attractor2,xyz0,tspan,p) ## problem is a Julia Object created
prob3 = ODEProblem(Lorenz_attractor3,xyz0,tspan,p)
# by the function ODEProblem - this object formulates the ODE proble

solution2 = solve(prob2)
solution3 = solve(prob3)  ## function solve, solves the ODE problem
### solution is a Julia Object


Lorenz_attractor1([0, 0, 0. ], xyz0, p, 10.)

Lorenz_attractor2([0, 0, 0. ], xyz0, p, 10.)


Lorenz_attractor3([0, 0, 0. ], xyz0, p, 10.)

plot(solution2,vars=(1,2,3)) ## variable 0 is time t, variable 1, 2 & 3 are X, Y Z, here we plot
plot(solution3,vars=(1,2,3)) ## variable 0 is time t, variable 1, 2 & 3 are X, Y Z, here we plot



### Example ODE of course

function odeEx1(dx,x,p,t)
        dx[1] = exp(t) - x[1] * sin(t)

end

tspan = (0., 10.)  # timespan
x0 = [2.]   # initial value
p =0

prob = ODEProblem(odeEx1,x0,tspan,p)
sol = solve(prob)

plot(sol)

b = 0.00018
β1 = 0.00414
β2 = 0.0115
mu = 4.563 * 10^-5
S = b/mu
P = 0
α1 = 0.10
α2 = 0.10
IA = 0
IS = 0
ψ = 0.0051
E = 0



function simulation(dS)
        dS = [b - (β1*S*P/ 1+α1 *P) - β2* S(IA +IS)/1+α2*(IA +IS) + ψE− muS]

end

dS = 0
tspan = (0.0,50.0)

prob4 = ODEProblem(simulation,dS,tspan)
