
function ApplyPointLoads(F, information, t = 0)
    """
        Simple code for the application of point loads to the load vector.
        DOCSTRING

    # Arguments:
    - `F`: Load vector. 
    - `information`: Struct containting all information for model.
    - `t`: Time [s], default is 0 seconds for cases that are not time dependent.
    """
    DOF = information.DOF
    for POINT_LOAD in information.BOUNDARY_CONDITIONS.POINT_LOADS
        NODE_ID = POINT_LOAD.NODE_ID
        DIRECTION = POINT_LOAD.DIRECTION
        F[(NODE_ID*DOF)+DIRECTION-1] += POINT_LOAD.MAGNITUDE(t)
    end

    return F
end


function Elimination(K, M, F, information, t = 0)
    """
        Application of elimination method to stifness matrix and load vector given a set 
        of single point constrains. 

    DOCSTRING

    # Arguments:
    - `K`: Stiffness matrix.
    - `M`: Mass matrix. 
    - `F`: Load vector. 
    - `information`: Struct containting all information for model.
    - `t`: Time [s].
    """

    BCS = information.BOUNDARY_CONDITIONS

    for DISPLACEMENT_CONDITION in BCS.DISPLACEMENT_CONDITIONS
        NODE_ID = DISPLACEMENT_CONDITION.NODE_ID
        DIRECTION = DISPLACEMENT_CONDITION.DIRECTION
        MAGNITUDE = DISPLACEMENT_CONDITION.MAGNITUDE
        OFFSET_VALUES = MAGNITUDE .* K[:, NODE_ID+DIRECTION-1]
        F -= OFFSET_VALUES
        K[(2*NODE_ID)-1+(DIRECTION-1), :] .= 0.0
        K[:, (2*NODE_ID)-1+(DIRECTION-1)] .= 0.0
        K[(2*NODE_ID)-1+(DIRECTION-1), (2*NODE_ID)-1+(DIRECTION-1)] = 1.0
        M[(2*NODE_ID)-1+(DIRECTION-1), :] .= 0.0
        M[:, (2*NODE_ID)-1+(DIRECTION-1)] .= 0.0
        M[(2*NODE_ID)-1+(DIRECTION-1), (2*NODE_ID)-1+(DIRECTION-1)] = 1.0
    end

    return K, M, F
end

function Penalty(K, M, F, information, c = 10e4, t = 0)
    """
        Application of penalty method to stifness matrix and load vector given a set of 
        single point constrains. Eventually this code will be updated to handle 
        multi-point constraints.

    DOCSTRING

    # Arguments:
    - `K`: Stiffness matrix.
    - `M`: Mass matrix. 
    - `F`: Load vector. 
    - `information`: Struct containting all information for model.
    - `c`: Method specific parameter. 
    - `t`: Time [s].
    """
    BCS = information.BOUNDARY_CONDITIONS
    c *= maximum(K)

    for DISPLACEMENT_CONDITION in BCS.DISPLACEMENT_CONDITIONS
        NODE_ID = DISPLACEMENT_CONDITION.NODE_ID
        DIRECTION = DISPLACEMENT_CONDITION.DIRECTION
        MAGNITUDE = DISPLACEMENT_CONDITION.MAGNITUDE
        K[(2*NODE_ID)-1+(DIRECTION-1), (2*NODE_ID)-1+(DIRECTION-1)] += c
        F[(2*NODE_ID)-1+(DIRECTION-1)] += c * MAGNITUDE
    end

    return K, M, F
end
