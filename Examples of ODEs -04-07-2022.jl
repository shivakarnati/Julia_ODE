using Plots # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)



function Ex1(dx, x, p, t)
    alpha = p[1]
    beta = p[2]
    dx[1] = alpha - beta *x[1]
    dx
end

tspan = (0., 50.)
x0 = [10.]
alpha = 1/3  ## every 3 days one individual
beta = 1/5  ## 5 days to recover
p = [alpha beta]

ode = ODEProblem(Ex1, x0, tspan, p)
sol = solve(ode)

plot(sol)



using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)



function Ex2(dx, x, p, t)
    gamma = p[1]
    lambda = p[2]
    dx[1] =  - gamma *x[1]
    dx[2] = gamma *x[1] - lambda * x[2]
    dx
end

tspan = (0., 50.)
x0 = [1000. 0.]
gamma = 1/1  ## every 1 day one individual
lambda = 1/7  ## 7 days to recover
p = [gamma lambda]

ode = ODEProblem(Ex2, x0, tspan, p)
sol = solve(ode)

plot(sol)





function Ex2(dx, x, p, t)
    dx[1] = -3*x[2]
    dx[2] = 3*x[1]
    dx
end

tspan = (0., 20.)
x0 = [3., 3.]
p = [ ]

ode1 = ODEProblem(Ex2, x0, tspan, p)
sol = solve(ode1)

plot(sol)


## diseases model


function Ex3(dx, x, p, t)
    dx[1] = -p[1]*x[1]*x[2] /p[3] + p[2]*x[2]
    dx[2] = p[1]*x[1]*x[2]/p[3] - p[2]*x[2]
    dx
end

beta = 0.5
gamma=1/7
N=10000
p = [beta gamma N]


tspan = (0., 20.)
x0 = [N-100, 100.]

ode1 = ODEProblem(Ex3, x0, tspan, p)
sol = solve(ode1)

plot(sol)


## predator prey model


alpha = 4 # growth rate of prey
beta = 0.1  # rate at which pedators feed on prey
delta= 0.1  # predation rate
gamma = 2 # death rate of predators

 function PP(dx, x, p, t)
    dx[1] = p[1]*x[1] - p[2]*x[1]*x[2]
    dx[2] = p[4]*x[1]*x[2] -p[3]*x[2]
    dx
end


p = [alpha, beta, gamma, delta]


tspan = (0., 100.)
x0 = [1000, 10.]

ode1 = ODEProblem(PP, x0, tspan, p)
sol = solve(ode1)

plot(sol)
