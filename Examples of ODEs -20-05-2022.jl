using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)



### Example 1
function LinODE(dx,x,p,t)
    dx[1] = -2.5 * x[1] + x[2]
    dx[2] = -0.5 * x[2]
end

tspan = (0.,20.)
x0 = [5., 53]
p = [ ]

ode = ODEProblem(LinODE,x0,tspan,p)
sol = solve(ode)

plot(sol)


###############

function ex2(dx, x, p, t)
    dx[1] = 3 * x[1] + x[2] - x[3]
    dx[2] = - x[2] + 3*x[3]
    dx[3] = -3 * x[3]
    dx
end
tspan = (0., 5.)
p = []
x0 = [5., 3., -5.]
prob = ODEProblem(ex2, x0, tspan, p)
sol = solve(prob)
plot(sol,ylims=[-10,10])

function ex3(dx, x, p, t)
    dx[1] = 1.5 * x[1] + 2*x[2]
    dx[2] = -1.3 * x[1] - x[2] + x[3]
    dx[3] = -x[3]
    dx
end
tspan = (0., 100.)
p = []
x0 = [1., 3., -1.]
prob = ODEProblem(ex3, x0, tspan, p)
sol = solve(prob)
plot(sol)


function ex4(dx, x, p, t)
    dx[1] = -1.5 * x[1] +  2* x[2]
    dx[2] = -1.3 * x[1] - x[2] + x[3]
    dx[3] = -x[3]
    dx
end
tspan = (0., 10.)
p = []
x0 = [1., 3., -1.]
prob = ODEProblem(ex4, x0, tspan, p)
sol = solve(prob)
plot(sol)



###############




### Example 2
function LinODE2(dx,x,p,t)
    dx[1] = 3 * x[1] + x[2] - x[3]
    dx[2] =  - x[2] + 3 * x[3]
    dx[3] = -3 * x[3]
end

tspan = (0.,20.)
x0 = [5., 3, -5]
p = [ ]

ode2 = ODEProblem(LinODE2,x0,tspan,p)
sol2 = solve(ode2)

plot(sol2)



### Example 3
function LinODE3(du,u,p,t)
    du[1] = 1.5 * u[1] + 2 * u[2]
    du[2] =  - 1.3* u[1] -  u[2] + u[3]
    du[3] = - u[3]
end

tspan = (0.,20.)
u0 = [1., 3, -1]
p = [ ]

ode3 = ODEProblem(LinODE3,u0,tspan,p)
sol3 = solve(ode3)

plot(sol3)



### Example 4
function LinODE4(du,u,p,t)
    du[1] = -1.5 * u[1] + 2 * u[2]
    du[2] =  - 1.3* u[1] -  u[2] + u[3]
    du[3] = - u[3]
end

tspan = (0.,20.)
u0 = [1., 3, -1]
p = [ ]

ode4 = ODEProblem(LinODE4,u0,tspan,p)
sol4 = solve(ode4)

plot(sol4)
