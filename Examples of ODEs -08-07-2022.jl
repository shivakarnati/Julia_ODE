using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)



function Ex2(dx, x, p, t)
    alpha, beta, gamma, mu = p
    dx[1] = alpha - (beta + gamma)  *x[1]   # dA/dt
    dx[2] = beta * x[1] - mu * x[2]  # dB/dt
    dx[3] = gamma * x[1] - mu * x[3]  # dC/dt
    dx
end

tspan = (0., 50.)
x0 = [10., 0., 0.]
alpha = 1/3  ## 3 days to fill up 1 unit
beta = 1/5  ## 5 days to be processed into unit B
gamma = 1/6  ## 6 to be processed into unit C
mu = 1/2  ## 2 days to process one unit

p = [alpha, beta, gamma, mu]

ode = ODEProblem(Ex2, x0, tspan, p)
sol = solve(ode)

plot(sol)

#### exercise 3

function Ex3(dx, x, p, t)
    lambda, gamma = p
    dx[1] =  - lambda * x[1]   # dS/dt
    dx[2] = lambda * x[1] - gamma * x[2]  # dI/dt
    dx[3] = gamma * x[2] # dR/dt
    dx
end

tspan = (0., 500.)
x0 = [10000., 0., 0.]
lambda = 1/120  ## 120 months until disease breaks out
gamma = 1/2  ## 2 months to recover

p = [lambda, gamma]

ode = ODEProblem(Ex3, x0, tspan, p)
sol = solve(ode)

plot(sol)

### SIR model version 2 - Force of infection


#### exercise 4

function Ex4(dx, x, p, t)
    beta, gamma, N = p
    dx[1] =  - beta/N * x[1] * x[2]   # dS/dt
    dx[2] = beta/N * x[1] * x[2] - gamma * x[2]  # dI/dt
    dx[3] = gamma * x[2] # dR/dt
    dx
end

N = 10000
tspan = (0., 12.)
x0 = [N-1, 1., 0.]
beta = 5  ## 5 contacts per months
gamma = 1/2  ## 2 months to recover

p = [beta, gamma, N]

ode = ODEProblem(Ex4, x0, tspan, p)
sol = solve(ode)


plot(sol)


### SIR model with Erlang states

function SIR_Erlang(dx, x, p, t)
    n = length(x) - 2  ### the number of Erlang states #compartments -2
    beta, gamma, N = p
    gamma = n * gamma ## pregression rate through Erlang
    Isum =  sum( x[2 : (n+1) ] )
    dx[1] =  - beta/N * x[1] * Isum  # dS/dt
    dx[2] = beta/N * x[1] * Isum  - gamma * x[2]  # dI1/dt
    for k = 3:(n+1)
        dx[k] = gamma * (x[k-1] - x[k]) # dI2/dt .... dIn/dt
    end
    dx[n+2] = gamma * x[n+1]
    dx
end


N = 10000
n = 15 # number of Erlang states
tspan = (0., 12.)
x0 = zeros(n+2)
x0[1] = N-1
x0[2] = 1
beta = 5  ## 5 contacts per months
gamma = 1/2  ## 2 months to recover

p = [beta, gamma, N]

ode = ODEProblem(SIR_Erlang, x0, tspan, p)
sol1 = solve(ode)


plot(sol1)


#### we want to plot the infections cummalitively
ttt = length(sol1.u)

out1=fill(0.,(ttt,3))
for t in 1 : ttt
    out1[t,1] = sol1.u[t][1]    ## susceptibless
    out1[t,2] = sum(sol1.u[t][2:(n+1)])   ## cmmulative infected
    out1[t,3] = sol1.u[t][n+2]  ## recovered
end

plot(plot(sol1.t,out1))

### 1 Erlang state - old model

N = 10000
n = 1 # number of Erlang states
tspan = (0., 12.)
x0 = zeros(n+2)
x0[1] = N-1
x0[2] = 1
beta = 5  ## 5 contacts per months
gamma = 1/2  ## 2 months to recover

p = [beta, gamma, N]

ode = ODEProblem(SIR_Erlang, x0, tspan, p)
sol = solve(ode)


plot(sol)


#### we want to plot the infections cummalitively
ttt = length(sol.u)

out=fill(0.,(ttt,3))
for t in 1 : ttt
    out[t,1] = sol.u[t][1]    ## susceptibless
    out[t,2] = sum(sol.u[t][2:(n+1)])   ## cmmulative infected
    out[t,3] = sol.u[t][n+2]  ## recovered
end

plot(plot(sol),plot(sol1.t,out1))



#### SIR model with birth and death


function SIR_BD(dx, x, p, t)
    N = sum(x)
    beta, gamma, mu = p
    dx[1] =  - beta/N * x[1] * x[2] + mu * ( N - x[1])   # dS/dt
    dx[2] = beta/N * x[1] * x[2] - ( gamma + mu ) * x[2]  # dI/dt
    dx[3] = gamma * x[2] - mu * x[3] # dR/dt
    dx
end

N = 10000
tspan = (0., 2200.)
x0 = [N-1, 1., 0.]
beta = 5  ## 5 contacts per months
gamma = 1/2  ## 2 months to recover
mu = 1/720   ## life expectancy is 60 years
p = [beta, gamma, mu]

ode = ODEProblem(SIR_BD, x0, tspan, p)
sol = solve(ode)


plot(sol)

## Equilibrium values
beta > gamma + mu


## endemic equilibrium

S=N*(gamma+mu)/beta
I=N*((1-(mu+gamma)/beta)/(1+gamma/mu))
R=N*((1-(mu+gamma)/beta)/(1+gamma/mu)) * gamma/mu
S+I+R
S

## First an epidemic outbreak and then euilibration. But fluctuations.
## Firs ta huge epidemic outbreak, then periodic epidemic outbreats that become less and less severe
## Eigenvalues of the endemic equilibrium are complex conjugate


## Change parameters


N = 10000
tspan = (0., 50.)
x0 = [N-1, 1., 0.]
beta = 5  ## 5 contacts per months
gamma = 1/2  ## 2 months to recover
mu = 1/12   ## life expectancy is 1 year
p = [beta, gamma, mu]

ode = ODEProblem(SIR_BD, x0, tspan, p)
sol = solve(ode)


plot(sol)

## now endemic equilibrium is reached rather fast


## Change parameters


N = 10000
tspan = (0., 500.)
x0 = [N-10, 10., 0.]
beta = 2  ## 5 contacts per months
gamma = 2  ## 1/2 months to recover
mu = 1/72   ## life expectancy is 72 months
p = [beta, gamma, mu]

ode = ODEProblem(SIR_BD, x0, tspan, p)
sol = solve(ode)


plot(sol)


### endemic outbreat and then we reach the disease free equilibrium
