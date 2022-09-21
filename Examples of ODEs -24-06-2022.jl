
#using Pkg
#Pkg.add("Polynomials")
using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)
using Polynomials
using LinearAlgebra

Polynomial([1,3,4,-1,2]) ## try the syntax

## polynomial of the for equilibria

pol = Polynomial([1,-21/20,-1/100,1/4000])
xx2 = roots(pol) ## these are the x2-components of equilibria xx^(2), xx^(3), xx^(4)

## x1- components of equilibria  xx^(2), xx^(3), xx^(4)
xx1= xx2 .* (1 .- xx2 ./ 20)

for k in 1:3
    x1 = xx1[k]
    x2 = xx2[k]
    dF =  [1- x1/5 -2*x2 ; 1  1- x2/10 ]
    println(eigvals(dF))
end

### All equilibria have at least one eigenvalue with pos. real part -> all are unstable




function Ex1(dx, x, p, t)
    dx[1] = x[1]*(1-x[1]/10)-x[2]^2
    dx[2] = x[2]*(1-x[2]/20)+x[1]
    dx
end


# First equilibrium xx^(1)=(0,0)
tspan = (0.,1)
#x0 = [xx1[1]+0.1, xx2[1]+0.1]
x0 = [0.1, 0.1]
p = [ ]

ode = ODEProblem(Ex1, x0, tspan, p)
sol = solve(ode)

plot(sol)
plot(sol,vars=(1,2))
plot(sol,vars=(1,2),xlims=(0.09,0.15),ylims=(0.09,0.15))



# Second  equilibrium xx^(2)
tspan = (0.,0.1)
x0 = [xx1[1]+0.1, xx2[1]+0.1]
#x0 = [0.1, 0.1]
p = [ ]

ode = ODEProblem(Ex1, x0, tspan, p)
sol = solve(ode)

plot(sol)
plot(sol,vars=(1,2))
#plot(sol,vars=(1,2),xlims=(0.09,0.15),ylims=(0.09,0.15))

## Solution diverges faster in x1 component, corresponding to the larger eigenvalues


# Third equilibrium xx^(3)
tspan = (0.,1.9)
x0 = [xx1[2]+0.01, xx2[2]-0.01]
#x0 = [0.1, 0.1]
p = [ ]

ode = ODEProblem(Ex1, x0, tspan, p)
sol = solve(ode)

plot(sol)
plot(sol,vars=(1,2))
#plot(sol,vars=(1,2),xlims=(0.09,0.15),ylims=(0.09,0.15))

## Solution diverges faster in x1 component, corresponding to the larger eigenvalues



function LogGro(dN, N, p, t)
    rho = p[1]
    K = p[2]
    dN[1] = rho * N[1] * (1- N[1]/K)
    dN
end


# First equilibrium xx^(1)=(0,0)
tspan = (0.,30)
#x0 = [xx1[1]+0.1, xx2[1]+0.1]
x0 = [1]
p = [0.5 , 10000 ]

ode = ODEProblem(LogGro, x0, tspan, p)
sol = solve(ode)

plot(sol)


x0 = [20000]
p = [0.5 , 10000 ]

ode = ODEProblem(LogGro, x0, tspan, p)
sol = solve(ode)

plot(sol)


### predator prey model with resource competition


function PPmodel(dx, x, p, t)
    lamx, lamy, Kx, Ky, gamma, beta = p
    dx[1] = x[1] * (lamx - lamx/Kx *x[1] - gamma* x[2])
    dx[2] = x[2] * (-lamy - lamy/Ky * x[2] + beta * x[1])
    dx
end



tspan = (0.,200)
#x0 = [xx1[1]+0.1, xx2[1]+0.1]
x0 = [100 100]
lamx = 10
lamy = 5
Kx =2000
Ky =1000
gamma =1.1
beta =0.1
p=[lamx, lamy, Kx, Ky, gamma, beta]
lamy/beta

(lamx^2/Kx+lamy*gamma)/((lamx*lamy)/(Kx*Ky)+gamma*beta)

ode = ODEProblem(PPmodel, x0, tspan, p)
sol = solve(ode)

plot(sol)
plot(sol,vars=(1,2))


## one arametrization without internal equilibrium


tspan = (0.,20)
#x0 = [xx1[1]+0.1, xx2[1]+0.1]
x0 = [100 200]
lamx = 10
lamy = 5
Kx =1000
Ky =1000
gamma =1.1
beta =0.004
p=[lamx, lamy, Kx, Ky, gamma, beta]
lamy/beta

(lamx^2/Kx+lamy*gamma)/((lamx*lamy)/(Kx*Ky)+gamma*beta)

ode = ODEProblem(PPmodel, x0, tspan, p)
sol = solve(ode)

plot(sol)
plot(sol,vars=(1,2))
