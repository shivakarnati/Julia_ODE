

using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)

### we need the function with 4 arguments du, u, t, and p

function ODE1(du,u,t,p)
    du = 0.8 * u
    du
end

tspan=(0.,10.)
u0 = [1.]
p=[]  ## we deine the parameter vector as empty

ODEprob = ODEProblem(ODE1,u0,tspan,p)

sol = solve(ODEprob)

plot(sol)


### different implementation

function ODE2(du,u,t,p)
    du = p * u
    du
end

tspan=(0.,10.)
u0 = [1.]
p=[0.8]  ## we deine the parameter vector as empty

ODEprob = ODEProblem(ODE2,u0,tspan,p)

sol = solve(ODEprob)

plot(sol)


### example SIR model

function SIR(dy,y,p,t)
    S, I, R = y
          # better might be dS, dI, dR = repeat([0], 3)
    dy[1] = -0.9 * S*I/(S+I+R)
    dy[2] =  0.9 * S*I/(S+I+R) - 1.2*I
    dy[3] = 1.2 * I
    dy
end

tspan = (0.,10.)
y0 = [10000., 100,0]
p = []
SIRsystem = ODEProblem(SIR,y0,tspan,p)
sol = solve(SIRsystem)

plot(sol)

## different implementation

function SIR1(dy,y,p,t)
    S, I, R = y
          # better might be dS, dI, dR = repeat([0], 3)

    dy[1] = -0.9 * S*I/(S+I+R)
    dy[2] =  0.9 * S*I/(S+I+R) - 1.2*I
    dy[3] = 1.2 * I

    dy
end

tspan = (0.,10.)
y0 = [10000., 100,0]
p = []
SIRsystem = ODEProblem(SIR1,y0,tspan,p)
sol = solve(SIRsystem)

plot(sol)


## another implementation
function SIR2(dy,y,p,t)
    S, I, R = y

    a1, b1 = p
    dS =  -a1 * S*I/(S+I+R)
    dI =  a1 * S*I/(S+I+R) - b1 *I
    dR = b1 * I
    dy .= [dS,dI,dR]
    dy
end

tspan = (0.,100.)
y0 = [10000., 1,0]
p = [2,  0.1]
SIRsystem1 = ODEProblem(SIR2,y0,tspan,p)
sol = solve(SIRsystem1)

plot(sol)

## for other parameters

tspan = (0.,100.)
y0 = [10000., 1,0]
p = [1.2,  0.3]
SIRsystem1 = ODEProblem(SIR2,y0,tspan,p)
sol = solve(SIRsystem1)

plot(sol)

## for other parameters

tspan = (0.,100.)
y0 = [10000., 1,0]
p = [0.8,  0.3]
SIRsystem1 = ODEProblem(SIR2,y0,tspan,p)
sol = solve(SIRsystem1)

plot(sol)

## for other parameters

tspan = (0.,150.)
y0 = [10000., 1,0]
p = [0.4,  0.3]    ## disease even less infectious
SIRsystem1 = ODEProblem(SIR2,y0,tspan,p)
sol = solve(SIRsystem1)

plot(sol)


### 2nd order ODE example rewritten as system of dim 2

function ODE3(dx,x,p,t)
    a = p
    dx[1] = x[2]
    dx[2] = 2 * x[2] - a * sin(x[1])
end

tspan=(0,5)
x0 = [1.,1]
p = 3.

ODEprob2 = ODEProblem(ODE3,x0,tspan,p)

sol2 = solve(ODEprob2)

plot(sol2)

## plot only solution if 2nd order ODE y(t)=x1


plot(sol2,vars=(0,1))
