using Plots
using DifferentialEquations


# model parameters
IA = 1000
IS = 1000
R = 2000
P = 5000
S = 2000
E = 3000


# parameters
    b = 1200000  #Birth rate of the human population
    mu = 0.4563  # Natural human death rate
    mu1 = 60   #Human life expectancy
    mu_p = 0.01724    #Natural death rate of pathogens in the environment
    #Life expectancy of pathogens in the environment
    alpha1 = 1.99    #Proportion of interaction with an infectious environment
    alpha2 = 1.99    #Proportion of interaction with an infectious individual
    beta1 = 0.00414  #Rate of transmission from S to E due to contact with P
    beta2 = 0.0115   #Rate of transmission from S to E due to contact with IA and/or IS
    delta = 0.7      #Proportion of symptomatic infectious people
    psi = 0.0051     #Progression rate from E back to S due to robust immune system
    w = 0.09         #Progression rate from E to either IA or IS
    sigma = 0.0018   #Death rate due to the coronavirus
    gamma_S = 0.05   #Rate of recovery of the symptomatic population
    gamma_A = 0.0714 #Rate of recovery of the asymptomatic human population
    eta1 = 0.1       #Rate of virus spread to environment by symptomatic infectious individuals
    eta2 = 0.05      #Rate of virus spread to environment by asymptomatic infectious individuals

x0 = [10]
alpha = 0.10
beta = 0.00414
p = [alpha1,alpha2,beta1,beta2,eta1,eta2,b,mu,mu_p]
psi = 0.0051

tspan = (0., 300)
p = [0, 0 , 0 ,0 , 0, 0, 0,0 ,0]

intial_Parameters = [IA,IS,R,P,E,S]

A = (beta1*S*P/(1+alpha1*P))
B = beta2*S*(IA+ IS)/(1+alpha2*(IA+IS))

function symptamatic(dxyz,initial_Parameters,p,tspan)

    IA, IS, R,P,E,S = initial_Parameters
    alpha1,alpha2,beta1,beta2,eta1,eta2,b,mu = p
    dxyz[1] = b - A - B + psi * E- mu*S
    print(dxyz[1])
    dxyz[2] = A + B - psi * E - mu*E - w*E

    dxyz[3] = (1-delta)*w*E - (mu + sigma)*IA - gamma1 * IA
    dxyz[4] = delta * w *E - (mu + sigma)*IS - mu * R
    dxyz[5] = gamma1*IS + gamma2*IA - mu*R
    dxyz[6] = eta1 * IA + eta2*IS - mu_p*P
    dxyz
end


 prob = ODEProblem(symptamatic,intial_Parameters,tspan,p)

solution = solve(prob)

plot(solution)
