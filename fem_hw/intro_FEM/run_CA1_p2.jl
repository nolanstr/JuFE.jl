include("../src/JuFE.jl")

information = JuFE.ReadYAML("./ca_1.yaml")
K, F = JuFE.GlobalAssembly(information)
K_PL, F_PL = JuFE.ApplyPointLoads(K, F, information)

K_Elim, F_Elim = JuFE.Elimination(copy(K_PL), copy(F_PL), information)
Q = JuFE.GaussElimination(K_Elim, F_Elim)
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)
internal_force = JuFE.T2D2_Internal_Force_From_Stress(information, stress)

println("USING ELIMINATION APPROACH:\n")
println("Displacement Under Load (Units --> m):\n")
for i in range(1, information.NX)
    print("Displacements (q_x, q_y) at Node $(i): $(reshape(Q[(2*i)-1:2*i], (1,2)))")
    println()
end

println("Stress:\n")
for i in range(1, 17)
    print(stress[i])
    println()
end

println("Total Force at Left Pin (fx, fy, units-->[N]): $(R[1:2])")

K_Pen, F_Pen = JuFE.Penalty(copy(K_PL), copy(F_PL), information)
Q = JuFE.GaussElimination(K_Pen, F_Pen)
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)
internal_force = JuFE.T2D2_Internal_Force_From_Stress(information, stress)

println("\nUSING PENALTY APPROACH:\n")
println("Displacement Under Load (Units --> m):\n")
for i in range(1, information.NX)
    print("Displacements (q_x, q_y) at Node $(i): $(reshape(Q[(2*i)-1:2*i], (1,2)))")
    println()
end

println("Stress:\n")
for i in range(1, 17)
    print(stress[i])
    println()
end

println("Total Force at Left Pin (fx, fy, units-->[N]): $(R[1:2])")

