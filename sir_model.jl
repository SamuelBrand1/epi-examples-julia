# Simple example of solving a simple SIR model and plotting it
using OrdinaryDiffEq, Plots

##

# Define the vector field function

function sir(du, u, p, t)
    # Get state
    S, I, R = u

    # Get the parameters
    R₀, γ = p

    # Set the dynamics
    du[1] = - γ * R₀ * S * I # dS
    du[2] =  γ * R₀ * S * I - γ * I # dI
    du[3] = γ * I # dR
end

f = ODEFunction(sir, syms = [:S, :I, :R])
# Wrap the vector field function in an ODEFunction. This can be used to add features, in this case I'm naming the variables.

# Define the initial conditions and parameters

p = [2.0, 0.1] # R0 is 2 and the mean infectious period is 10 days
u0 = [0.999, 0.001, 0.0] # Model is in proportion of population, this initialises with 0.1% of population infected.

## Use OrdinaryDiffEq to define and solve the ODE model between time 0 and day 150

prob = ODEProblem(f, u0, (0.0, 150.0), p)
sol = solve(prob, Tsit5())

## Use Plots to plot the solution
# `plot` dispatches on the input being an ODESolution object.
plot(sol)


