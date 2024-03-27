
include("global_assembly.jl")
include("apply_boundary_conditions.jl")

function NonLinearSolver(information, tf, dt)
    """
        Currently this code does not consider dampening!

    DOCSTRING

    # Arguments:
    - `information`: Struct containting all information for model.
    - `tf`: Final time of simulation [s].
    - `dt`: Time step [s].
    """
    if information.METHOD == "Central Difference"
        method = CentralDifference
    elseif information.METHOD == "Trapezoidal Rule"
        method = TrapezoidalRule
    elseif information.METHOD == "Wilson Theta"
        method = WilsonTheta
    elseif information.Method == "Newmark Beta"
        method = NewmarkBeta
    else
        println(
            "No valid Non-Linear Solver method provided, defaulting to Central Difference!",
        )
        method = CentralDifference
    end

    if information.BC_METHOD == "Elimination"
        BC_method = Elimination
    elseif information.BC_METHOD == "Penalty"
        BC_method = Penalty
    else
        println(
            "No valid Boundary Condition method provided, defaulting to Elimination Method!",
        )
        BC_method = Elimination
    end

    total_dofs = information.NX * information.DOF
    method_args = information.METHOD_ARGS

    time = 0:dt:tf
    steps = length(time)
    u = Matrix{Float64}(undef, total_dofs, steps)
    u[:, 1] = zeros(total_dofs)
    v = Matrix{Float64}(undef, total_dofs, steps)
    v[:, 1] = zeros(total_dofs)
    a = Matrix{Float64}(undef, total_dofs, steps)
    a[:, 1] = zeros(total_dofs)

    for (i, t) in enumerate(time[1:end-1])

        K, M, F = NonLinearGlobalAssembly(information, t)
        F = ApplyPointLoads(F, information, t)
        K, M, F = BC_method(K, M, F, information, t)
        u_tn, v_tn, a_tn = method(K, M, F, u, v, a, dt, i, method_args)
        u[:, i+1] = u_tn

    end

    return u, v, a
end


function CentralDifference(K, M, F, u, v, a, dt, i, method_args)
    """
        CentralDifference(K, M, F, u, v, a, dt, i, method_args)

    DOCSTRING

    # Arguments:
    - `K`: Stiffness matrix.
    - `M`: Mass matrix.
    - `F`: Load Vector.
    - `u`: Displacement vector(s).
    - `v`: Velocity vector(s).
    - `a`: Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`: Iteration.
    - `method_args`: Method specific optional arguments.
    """
    u_curr = u[:, i]
    if i == 1
        u_prev = u[:, 1] * 0 #Currently just setting to zero before I figure a better way to deal with this.
    else
        u_prev = u[:, i-1]
    end
    left_M_term = M / (dt^2)
    u_curr_terms = (K - ((2 / (dt^2)) * M)) * u_curr
    u_prev_terms = (M / (dt^2)) * u_prev
    effective_load = F - u_curr_terms - u_prev_terms
    u_tn = inv(left_M_term) * effective_load
    v_tn = (u_tn - u_curr) / dt
    a_tn = (v_tn)
    return u_tn, v_tn, a_tn
end

function TrapezoidalRule(K, M, F, u, v, a, dt, i, method_args)
    """
        TrapezoidalRule(K, M, F, u, v, a, dt, i, method_args)

    DOCSTRING

    # Arguments:
    - `K`: Stiffness matrix.
    - `M`: Mass matrix.
    - `F`: Load Vector.
    - `u`: Displacement vector(s).
    - `v`: Velocity vector(s).
    - `a`: Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`: Iteration.
    - `method_args`: Method specific optional arguments.
    """
    return u_tn, v_tn, a_tn
end

function WilsonTheta(K, M, F, u, v, a, dt, i, method_args)
    """
        WilsonTheta(K, M, F, u, v, a, dt, i, method_args)

    DOCSTRING

    # Arguments:
    - `K`: Stiffness matrix.
    - `M`: Mass matrix.
    - `F`: Load Vector.
    - `u`: Displacement vector(s).
    - `v`: Velocity vector(s).
    - `a`: Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`: Iteration.
    - `method_args`: Method specific optional arguments.
    """
    return u_tn, v_tn, a_tn
end

function NewmarkBeta(K, M, F, u, v, a, dt, i, method_args)
    """
        NewmarkBeta(K, M, F, u, v, a, dt, i, method_args)

    DOCSTRING

    # Arguments:
    - `K`: Stiffness matrix.
    - `M`: Mass matrix.
    - `F`: Load Vector.
    - `u`: Displacement vector(s).
    - `v`: Velocity vector(s).
    - `a`: Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`: Iteration.
    - `method_args`: Method specific optional arguments.
    """
    alpha, beta = method_args
    u_curr = u[:, i]
    if i == 1
        u_prev = u[:, 1] * 0
    else
        u_prev = u[:, i-1]
    end
    left_M_term = M / (dt^2)
    u_curr_terms = (K - ((2 / (dt^2)) * M)) * u_curr
    u_prev_terms = (M / (dt^2)) * u_prev
    effective_load = F - u_curr_terms - u_prev_terms
    u_tn = inv(left_M_term) * effective_load
    v_tn = (u_tn - u_curr) / dt
    a_tn = (v_tn)
    return u_tn, v_tn, a_tn
end
