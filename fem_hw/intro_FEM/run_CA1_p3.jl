include("../src/JuFE.jl")

F3_7 = 7071.
F7_9 = -7071.

F5_10 = -11180.0

analytical_force = [F3_7 F7_9 F5_10]
analytical_force = reshape(analytical_force, (:,1))

information = JuFE.ReadYAML("./ca_1.yaml")
K, F = JuFE.GlobalAssembly(information)
K_PL, F_PL = JuFE.ApplyPointLoads(K, F, information)

#(a) Elimination Method
K_Elim, F_Elim = JuFE.Elimination(copy(K_PL), copy(F_PL), information)
Q = JuFE.GaussElimination(K_Elim, F_Elim)
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)
internal_force = JuFE.T2D2_Internal_Force_From_Stress(information, stress)
check_terms = reshape([internal_force[9, 1],
               internal_force[14, 1],
               internal_force[11, 1]], (:,1))
INDIVIDUAL_ERROR = broadcast(abs,
                (analytical_force-check_terms) ./ analytical_force)
TOTAL_ERROR = abs((sum(analytical_force)-sum(check_terms))/sum(analytical_force))
println("Elimination Method:")
println("Relative Errors in right side non-horizontal members: $(INDIVIDUAL_ERROR)")
println("Relative Total Error in right side non-horizontal members: $(TOTAL_ERROR)")
#display(internal_force)

#(b) Penalty method with C = max(K)*10000
K_Pen1, F_Pen1 = JuFE.Penalty(copy(K_PL), copy(F_PL), information, 10e4)
Q = JuFE.GaussElimination(K_Pen1, F_Pen1)
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)
internal_force = JuFE.T2D2_Internal_Force_From_Stress(information, stress)

check_terms = reshape([internal_force[9, 1],
               internal_force[14, 1],
               internal_force[11, 1]], (:,1))
INDIVIDUAL_ERROR = broadcast(abs,
                (analytical_force-check_terms) ./ analytical_force)
TOTAL_ERROR = abs((sum(analytical_force)-sum(check_terms))/sum(analytical_force))
println("\nPenalty method with C=10000*K_max:")
println("Relative Errors in right side non-horizontal members: $(INDIVIDUAL_ERROR)")
println("Relative Total Error in right side non-horizontal members: $(TOTAL_ERROR)")
#
#(c) Penalty method with C = max(K)*100
K_Pen2, F_Pen2 = JuFE.Penalty(copy(K_PL), copy(F_PL), information, 10e2)
Q = JuFE.GaussElimination(K_Pen2, F_Pen2)
stress = JuFE.T2D2_Stress(information, Q)
R = JuFE.Reaction(K, Q, F)
internal_force = JuFE.T2D2_Internal_Force_From_Stress(information, stress)

check_terms = reshape([internal_force[9, 1],
               internal_force[14, 1],
               internal_force[11, 1]], (:,1))
INDIVIDUAL_ERROR = broadcast(abs,
                (analytical_force-check_terms) ./ analytical_force)
TOTAL_ERROR = abs((sum(analytical_force)-sum(check_terms))/sum(analytical_force))
println("\nPenalty method with C=100*K_max:")
println("Relative Errors in right side non-horizontal members: $(INDIVIDUAL_ERROR)")
println("Relative Total Error in right side non-horizontal members: $(TOTAL_ERROR)")

