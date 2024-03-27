using LinearAlgebra

function T2D2_Stress(information, displacements)
    """
    Code to compute the T2D2 element stress. This code currently performs this over
    all elements in information but eventually this will need to be made into a 
    function that performs this for a single element and is called by a global
    stress analysis function that can take a mixture of element types.
    """
    stress = zeros(information.NELX)
    for element in information.ELEMENTS
        E = information.MATERIALS[element.MATERIAL_ID].E
        ID = element.ID
        NODE1, NODE2 = element.NODE_IDS[1], element.NODE_IDS[2]
        IDXS = [2 * NODE1 - 1, 2 * NODE1, 2 * NODE2 - 1, 2 * NODE2]
        NODE1 = information.NODAL_COORDS[NODE1]
        NODE2 = information.NODAL_COORDS[NODE2]
        le = sqrt((NODE2.x - NODE1.x)^2 + (NODE2.y - NODE1.y)^2)
        l = (NODE2.x - NODE1.x) / le
        m = (NODE2.y - NODE1.y) / le
        q = reshape(displacements[IDXS], (1, 4))
        stress[ID] = (E / le) * dot(Matrix([-l -m l m]), q)
    end
    return reshape(stress, (:, 1))
end

function T2D2_Strain(information, displacements)
    """
    Code to compute the T2D2 element strain. This code currently performs this over
    all elements in information but eventually this will need to be made into a 
    function that performs this for a single element and is called by a global
    strain analysis function that can take a mixture of element types.
    """
    strain = zeros(information.NELX)
    for element in information.ELEMENTS
        E = information.MATERIALS[element.MATERIAL_ID].E
        ID = element.ID
        NODE1, NODE2 = element.NODE_IDS[1], element.NODE_IDS[2]
        IDXS = [2 * NODE1 - 1, 2 * NODE1, 2 * NODE2 - 1, 2 * NODE2]
        NODE1 = information.NODAL_COORDS[NODE1]
        NODE2 = information.NODAL_COORDS[NODE2]
        le = sqrt((NODE2.x - NODE1.x)^2 + (NODE2.y - NODE1.y)^2)
        l = (NODE2.x - NODE1.x) / le
        m = (NODE2.y - NODE1.y) / le
        q = reshape(displacements[IDXS], (1, 4))
        strain[ID] = dot(Matrix([-l -m l m]), q)
    end
    return reshape(strain, (:, 1))
end

function T2D2_Internal_Force_From_Stress(information, stress)
    """
    Computes the internal force of a T2D2 element from the stress.
    """
    force = zeros(information.NELX)
    for element in information.ELEMENTS
        E = information.MATERIALS[element.MATERIAL_ID].E
        A = element.AREA
        ID = element.ID
        force[ID] = stress[ID] * A
    end
    return reshape(force, (:, 1))
end

function T2D2_Internal_Force_From_Strain(information, strain)
    """
    Computes the internal force of a T2D2 element from the strain.
    """
    force = zeros(information.NELX)
    for element in information.ELEMENTS
        E = information.MATERIALS[element.MATERIAL_ID].E
        A = element.AREA
        ID = element.ID
        force[ID] = (E * strain[ID]) * A
    end
    return reshape(force, (:, 1))
end

function Reaction(K, Q, F)
    """
    Computes the reaction force at each node in the FE model.
    """
    n, _ = size(K)
    R = (K * reshape(Q, (n, 1))) - reshape(F, (n, 1))
    return R
end
