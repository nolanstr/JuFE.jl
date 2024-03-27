using Plots
include("../../src/JuFE.jl")

"""
Time parameters used for analytica and Central Difference.
"""
tf = 0.01 
dt =1e-3
"""
Analytical Solution Using Hooke's Law in 1D.
"""
LOAD = t -> (0.5*t + 0.3*t^2) * 10^2
time = 0:dt:tf
L = 1
E = 10^6
A = 1
u_analytical = (L/(E*A)) * map(LOAD, time)
p = plot(time, u_analytical, label="Analytical")
xlabel!("Time [s]")
ylabel!("Displacement [m]")


information = JuFE.ReadYAML("./hw1_p4.yaml")
u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(time, u[2,:], label="Central Difference", ls=:dot)
savefig(p, "myplot.pdf")
