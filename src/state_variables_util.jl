
function initializeStateVariables(information, BC_method, steps)
    ndofs = information.NX * information.DOF
    K_0, M_0, F_0 = NonLinearGlobalAssembly(information, 0)
    F_0 = ApplyPointLoads(F_0, information, 0)
    K_0, M_0, F_0 = BC_method(K_0, M_0, F_0, information, 0)
    K = Array{Float64}(undef, steps, ndofs, ndofs)
    M = Array{Float64}(undef, steps, ndofs, ndofs)
    F = Matrix{Float64}(undef, steps, ndofs)
    K[1, :, :] = K_0
    M[1, :, :] = M_0
    F[1, :] = F_0
    return K, M, F
end

function updateStateVariables(K, M, F, information, t, i, BC_method)
    K_i, M_i, F_i = NonLinearGlobalAssembly(information, t)
    F_i = ApplyPointLoads(F_i, information, t)
    K_i, M_i, F_i = BC_method(K_i, M_i, F_i, information, t)
    K[i, :, :] = K_i
    M[i, :, :] = M_i
    F[i, :] = F_i
    return K, M, F
end

function initializeSolutionFields(information, steps)
    ndofs = information.NX * information.DOF
    u = Matrix{Float64}(undef, steps, ndofs)
    v = Matrix{Float64}(undef, steps, ndofs)
    a = Matrix{Float64}(undef, steps, ndofs)
    u[1, :] = zeros(ndofs)
    v[1, :] = zeros(ndofs)
    a[1, :] = zeros(ndofs)
    return u, v, a
end

function updateSolutionFields(K, M, F, u, v, a, dt, i, method, method_args)
    u_tn, v_tn, a_tn = method(K, M, F, u, v, a, dt, i, method_args)
    u[i, :] = u_tn
    v[i, :] = v_tn
    a[i, :] = a_tn
    return u, v, a
end
