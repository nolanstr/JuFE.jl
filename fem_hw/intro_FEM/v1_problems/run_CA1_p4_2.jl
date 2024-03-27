include("../src/JuFE.jl")

information = JuFE.ReadYAML("./ca_1_p4_2.yaml")
K, F = JuFE.GlobalAssembly(information)
K_PL, F_PL = JuFE.ApplyPointLoads(K, F, information)

K_Elim, F_Elim = JuFE.Elimination(copy(K_PL), copy(F_PL), information)
Q = JuFE.GaussElimination(K_Elim, F_Elim)
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)
internal_force = JuFE.T2D2_Internal_Force_From_Stress(information, stress)
println("Internal Forces:")
display(internal_force)

println("\nUSING ELIMINATION APPROACH:")
println("Total Force at Left Pin: $(R[1:2])")

K_Pen, F_Pen = JuFE.Penalty(copy(K_PL), copy(F_PL), information)
Q = JuFE.GaussElimination(K_Pen, F_Pen)
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)
internal_force = JuFE.T2D2_Internal_Force_From_Stress(information, stress)

println("\nUSING PENALTY APPROACH:")
println("Total Force at Left Pin: $(R[1:2])")
