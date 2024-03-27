
module MassModule

function T2D2_CONSISTENT_MASS(MATERIAL, AREA, NODE1, NODE2)
    """
        T2D2_CONSISTENT_MASS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """
    M = Matrix([[0 0 0 0]; [0 0 0 0]])
    return M
end

function T2D2_HRZ_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)
    """
        T2D2_HRZ_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """
    M = Matrix([[0 0 0 0]; [0 0 0 0]])
    return M
end

function T2D2_OPTIMAL_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)
    """
        T2D2_OPTIMAL_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """
    M = Matrix([[0 0 0 0]; [0 0 0 0]])
    return M
end

function B11_CONSISTENT_MASS(MATERIAL, AREA, NODE1, NODE2)
    """
        B11_CONSISTENT_MASS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """
    le = NODE2.x - NODE1.x
    M = (MATERIAL.rho * AREA * le / 2) * Matrix([[2 1]; [1 2]])
    return M
end

function B11_HRZ_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)
    """
        B11_HRZ_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """
    le = NODE2.x - NODE1.x
    M = (MATERIAL.rho * AREA * le / 2) * Matrix([[1 0]; [0 1]])
    return M
end

function B11_OPTIMAL_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)
    """
        B11_OPTIMAL_LUMPED_MASS(MATERIAL, AREA, NODE1, NODE2)

    DOCSTRING

    # Arguments:
    - `MATERIAL`: Material dictionary containing material properties. 
    - `AREA`: Cross-sectional area of element.
    - `NODE1`: Node dictionary of element local node number 1.
    - `NODE2`: Node dictionary of element local node number 2.
    """
    le = NODE2.x - NODE1.x
    sum_evals = Matrix([[1 0]; [0 1]])
    M = (MATERIAL.rho * AREA * le / 2) * sum_evals
    return M
end

end
