include("global_assembly.jl")
include("apply_boundary_conditions.jl")
include("state_variables_util.jl")

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
        method_args = information.METHOD_ARGS
    elseif information.METHOD == "Trapezoidal Rule"
        method = TrapezoidalRule
        method_args = information.METHOD_ARGS
    elseif information.METHOD == "Wilson Theta"
        method = WilsonTheta
        method_args = [1.4]
    elseif information.METHOD == "Newmark Beta"
        method = NewmarkBeta
        method_args = information.METHOD_ARGS
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

    time = 0:dt:tf
    steps = length(time)

    K, M, F = initializeStateVariables(information, BC_method, steps)
    u, v, a = initializeSolutionFields(information, steps)

    for i in range(2, length(time))
        t = time[i]
        K, M, F = updateStateVariables(K, M, F, information, t, i, BC_method)
        u, v, a = updateSolutionFields(K, M, F, u, v, a, dt, i, method, method_args)
    end

    return time, u, v, a
end


function CentralDifference(K, M, F, u, v, a, dt, i, method_args)
    """
        CentralDifference(K, M, F, u, v, a, dt, i, method_args)
        This function acts as a wrapper for Newmark-Beta with
        parameters beta = 0.0 and gamma = 0.5.
    DOCSTRING

    # Arguments:
    - `K`:  Stiffness matrix.
    - `M`:  Mass matrix.
    - `F`:  Load Vector.
    - `u`:  Displacement vector(s).
    - `v`:  Velocity vector(s).
    - `a`:  Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`:  Iteration.
    - `method_args`: Method specific optional arguments.
    """
    return NewmarkBeta(K, M, F, u, v, a, dt, i, [1E-10, 0.5])
end

function TrapezoidalRule(K, M, F, u, v, a, dt, i, method_args)
    """
        TrapezoidalRule(K, M, F, u, v, a, dt, i, method_args)
        This function acts as a wrapper for Newmark-Beta with
        parameters beta = 0.25 and gamma = 0.5.

    DOCSTRING

    # Arguments:
    - `K`:  Stiffness matrix.
    - `M`:  Mass matrix.
    - `F`:  Load Vector.
    - `u`:  Displacement vector(s).
    - `v`:  Velocity vector(s).
    - `a`:  Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`:  Iteration.
    - `method_args`: Method specific optional arguments.
    """
    return NewmarkBeta(K, M, F, u, v, a, dt, i, [0.25, 0.5])
end

function WilsonTheta(K, M, F, u, v, a, dt, i, method_args)
    """
        WilsonTheta(K, M, F, u, v, a, dt, i, method_args)

    DOCSTRING

    # Arguments:
    - `K`:  Stiffness matrix.
    - `M`:  Mass matrix.
    - `F`:  Load Vector.
    - `u`:  Displacement vector(s).
    - `v`:  Velocity vector(s).
    - `a`:  Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`:  Iteration.
    - `method_args`: Theta parameter specific to Wilson-Theta method.
    """
    theta = method_args[1]

    a_0 = 6 / (theta * dt)^2
    a_1 = 3 / (theta * dt)
    a_2 = 2 * a_1
    a_3 = (theta * dt) / 2
    a_4 = a_0 / theta
    a_5 = -a_2 / theta
    a_6 = 1 - (3 / theta)
    a_7 = dt / 2
    a_8 = (dt^2) / 6

    u_t = u[i-1, :]
    v_t = v[i-1, :]
    a_t = a[i-1, :]

    effective_stiffness = K[i, :, :] + a_0 * M[i, :, :] #+ a_1 * C
    effective_load = (
        F[i-1, :] +
        theta * (F[i, :] - F[i-1, :]) +
        M[i, :, :] * (a_0 * u_t + a_2 * v_t + 2 * a_t)
    )
    #+ C * (a_1 * u_t + 2 * v_t + a_3 * a_t) 
    u_t_theta = effective_stiffness \ effective_load

    a_tn = a_4 * (u_t_theta - u_t) + a_5 * v_t + a_6 * a_t
    v_tn = v_t + a_7 * (a_tn + a_t)
    u_tn = u_t + dt * v_t + a_8 * (a_tn + 2 * a_t)

    return u_tn, v_tn, a_tn
end

function NewmarkBeta(K, M, F, u, v, a, dt, i, method_args)
    """
        NewmarkBeta(K, M, F, u, v, a, dt, i, method_args)

    DOCSTRING

    # Arguments:
    - `K`:  Stiffness matrix.
    - `M`:  Mass matrix.
    - `F`:  Load Vector.
    - `u`:  Displacement vector(s).
    - `v`:  Velocity vector(s).
    - `a`:  Acceleration vector(s).
    - `dt`: Time step [s].
    - `i`:  Iteration.
    - `method_args`: Alpha and beta parameters for Newmark-Beta method. 
    """
    alpha = method_args[1]
    beta = method_args[2]

    a_0 = 1 / (alpha * (dt^2))
    a_1 = beta / (alpha * dt)
    a_2 = 1 / (alpha * dt)
    a_3 = (1 / (2 * alpha)) - 1
    a_4 = (beta / alpha) - 1
    a_5 = (dt / 2) * ((beta / alpha) - 2)
    a_6 = dt * (1 - beta)
    a_7 = beta * dt

    u_t = u[i-1, :]
    v_t = v[i-1, :]
    a_t = a[i-1, :]

    effective_stiffness = K[i, :, :] + (a_0 * M[i, :, :]) #+ a_1 * C[i]
    effective_load = F[i, :] + (M[i, :, :] * ((a_0 * u_t) + (a_2 * v_t) + (a_3 * a_t)))
    #+ C[i] * (a_1 * u_t + a_4 * v_t + a_5 * a_t) 
    u_tn = effective_stiffness \ effective_load
    a_tn = (a_0 * (u_tn - u_t)) - (a_2 * v_t) - (a_3 * a_t)
    v_tn = v_t + (a_6 * a_t) + (a_7 * a_tn)

    return u_tn, v_tn, a_tn
end
