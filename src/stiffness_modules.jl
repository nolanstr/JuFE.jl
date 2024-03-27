
module StiffnessModule

function T2D2_STIFFNESS(MATERIAL, AREA, NODE1, NODE2)
    """
        T2D2_STIFFNESS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """
    le = sqrt((NODE2.x - NODE1.x)^2 + (NODE2.y - NODE1.y)^2)
    l = (NODE2.x - NODE1.x) / le
    m = (NODE2.y - NODE1.y) / le
    k = ((MATERIAL.E * AREA) / le) .* Matrix([[1.0 -1.0]; [-1.0 1.0]])
    L = Matrix([[l m 0 0]; [0 0 l m]])
    k = L' * k * L
    f = L' * [0; 0] #Assuming no body force for this case.
    return le, k, f
end

function B11_STIFFNESS(MATERIAL, AREA, NODE1, NODE2)
    """
        B11_STIFFNESS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """

    le = NODE2.x - NODE1.x
    k = ((MATERIAL.E * AREA) / le) .* Matrix([[1.0 -1.0]; [-1.0 1.0]])
    f = k * [0; 0] #Assuming no body force for this case.
    return le, k, f
end

end
