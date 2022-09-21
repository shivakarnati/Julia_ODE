
using Pkg
#Pkg.add("Polynomials")
using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)
using Polynomials
using LinearAlgebra

Polynomial([1,3,4,-1,2]) ## try the syntax

## polynomial of the for equilibria

pol = Polynomial([1,-21/20,-1/100,1/4000])
xx2 = roots(pol)


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

tspan = (0.,0.02)
x0 = [xx1[1], xx2[1]+0.1]
p = [ ]

ode = ODEProblem(Ex1, x0, tspan, p)
sol = solve(ode)

plot(sol)
