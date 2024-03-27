include("../../src/JuFE.jl")

information = JuFE.ReadYAML("./ca_1.yaml")
K, F = JuFE.GlobalAssembly(information)
K_PL, F_PL = JuFE.ApplyPointLoads(K, F, information)

K_Elim, F_Elim = JuFE.Elimination(copy(K_PL), copy(F_PL), information)
println(K-K_Elim)
Q  = K_Elim \ F_Elim
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)

print("Displacements:")
display(Q)
print("Stress:")
display(stress)
print("Reaction:")
display(R)
