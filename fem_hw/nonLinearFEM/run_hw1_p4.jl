
using Plots
using LaTeXStrings
include("../../src/JuFE.jl")
YAML_DIR = "./YAML_FILES/"

"""
Time parameters used for analytica and Central Difference.
"""
tf =0.01
dt = 1e-3
dt_consistent = 0.0005484827557301444 
dt_lumped = 0.00095
"""
Analytical Solution Using Hooke's Law in 1D.
"""
LOAD = t -> (0.5 * t + 0.3 * t^2) * 10^2
time = 0:dt:tf
L = 1
E = 10^6
A = 1
u_analytical = (L / (E * A)) * map(LOAD, time)
p = plot(time, u_analytical, label = "Analytical", color = :black)
xlabel!("Time [s]")
ylabel!("Displacement [m]")

STEPS = 10

# Central Difference
information = JuFE.ReadYAML(YAML_DIR * "/hw1_p4_central_difference_consistent.yaml")
dt = dt_consistent
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Central Difference (Consistent Mass Matrix)",
    ls = :dot,
    color = :red,
)

information = JuFE.ReadYAML(YAML_DIR * "hw1_p4_central_difference_lumped.yaml")
dt = dt_lumped
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Central Difference (HRZ Lumped Mass Matrix)",
    ls = :dashdot,
    color = :red,
)

# Trapezoidal Rule 
information = JuFE.ReadYAML(YAML_DIR * "hw1_p4_trapezoidal_rule_consistent.yaml")
dt = dt_consistent 
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Trapezoidal Rule (Consistent Mass Matrix)",
    ls = :dot,
    color = :blue,
)

information = JuFE.ReadYAML(YAML_DIR * "hw1_p4_trapezoidal_rule_lumped.yaml")
dt = dt_lumped
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Trapezoidal Rule (HRZ Lumped Mass Matrix)",
    ls = :dashdot,
    color = :blue,
)

# Wilson-Thete 
information = JuFE.ReadYAML(YAML_DIR * "hw1_p4_wilson_theta_consistent.yaml")
dt = dt_consistent 
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Wilson-" * L"\Theta" * " (Consistent Mass Matrix)",
    ls = :dot,
    color = :green,
)

information = JuFE.ReadYAML(YAML_DIR * "hw1_p4_wilson_theta_lumped.yaml")
dt = dt_lumped
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Wilson-" * L"\Theta" * " (HRZ Lumped Mass Matrix)",
    ls = :dashdot,
    color = :green,
)

# Newmark-Beta
information = JuFE.ReadYAML(YAML_DIR * "hw1_p4_newmark_beta_consistent.yaml")
dt = dt_consistent 
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Newmark-" * L"\beta" * " (Consistent Mass Matrix)",
    ls = :dot,
    color = :orange,
)

#information = JuFE.ReadYAML(YAML_DIR * "hw1_p4_newmark_beta_lumped.yaml")
dt = dt_lumped 
#tf = dt*STEPS
time, u, v, a = JuFE.NonLinearSolver(information, tf, dt)
plot!(
    time,
    u[:, 2],
    label = "Newmark-" * L"\beta" * " (HRZ Lumped Mass Matrix)",
    ls = :dashdot,
    color = :orange,
)

savefig(p, "time_displacement_plots.pdf")
