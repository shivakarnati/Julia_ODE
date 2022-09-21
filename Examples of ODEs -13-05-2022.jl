using Plots  # load package plots (Documentation: https://docs.juliaplots.org/stable/)
using DifferentialEquations # (Documentation: https://diffeq.sciml.ai/stable/)
                            # (Further: https://juliapackages.com/p/differentialequations)



## exact solution
function LogGrowth(t,K,u0,gamma)
    U0 = u0/(u0-K)
    U0 .*K .* exp.(gamma .*t) ./(U0 .* exp.(gamma .*t) .-1)
end

## Assume time in years
# birth rate = 2 individuals/ year /individual
# death rate = 1 individual /year /individual
# u0 ... initial pop size = 100 individuals

## number of individuals after 10 years almost 100 000
LogGrowth(10, 100000,100,1)
## number of individuals after 20 years = 100 000
LogGrowth(20, 100000,100,1)




plot((1:2000) ./100,LogGrowth((1:2000) ./100, 100000,100,1))

#### optimize the code - we want to programm efficiently, so save on multplications


## exact solution
function LogGrowth1(t,K,u0,gamma)
    #U0 .*K .* exp.(gamma .*t) ./(U0 .* exp.(gamma .*t) .-1)
    K  ./(1 .- exp.( - gamma .*t) .* (u0-K) ./u0 )
end

tt =(1:2000) ./100

plot(tt,LogGrowth1(tt, 100000,100,1))

### What happens if the pop exceeds carrying capacity K initially?


plot(tt,[LogGrowth1(tt, 100000,200000,1),LogGrowth1(tt, 100000,200,1)])


### Solving the ODE numericallly


### make sure to deine the operations element wise! Use . Operator!

function LGODE(du,u,p,t)
    K, gamma = p
    du .= gamma .* u .* (1 .- u/K)
    du
end

tspan = (0.,20.)
u0 = [200.]
gamma = 1.
K = 100000.
p = [K, gamma]

ode = ODEProblem(LGODE,u0,tspan,p)
sol = solve(ode)

plot(sol)
