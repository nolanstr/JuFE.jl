
function GlobalAssembly(information)
    """
        Simple code for assembling global stifness and load vector given T2DE elements.
        This code will be updated to handle different elements and a micture of elements.
        Approach: Generate list of lists where each index in the list relates to the i^th
        element and obtains the idxs of that elements nodes and use these to perform the 
        _i/_j loop.

    DOCSTRING

    # Arguments:
    - `information`: Struct containting all information for model.
    """
    K = zeros(information.NX * 2, information.NX * 2)
    F = zeros(information.NX * 2, 1)

    for element in information.ELEMENTS

        #hard coding this section for now for 4x4 element stiffness
        i, j = element.NODE_IDS[1], element.NODE_IDS[2]
        idxs = [2 * i - 1, 2 * i, 2 * j - 1, 2 * j]
        for _i in range(1, 4)
            for _j in range(1, 4)
                K[idxs[_i], idxs[_j]] += element.K[_i, _j]
            end
        end

        F[2*i-1] += element.LOAD[1]
        F[2*i] += element.LOAD[2]
        F[2*j-1] += element.LOAD[3]
        F[2*j] += element.LOAD[4]
    end
    return (K, F)
end

function NonLinearGlobalAssembly(information, t)
    """
        This is a copy of GlobalAssembly but will be further updated in the future 
        to capture nonlinearities that arrise in the stiffness and mass matrices. 

    DOCSTRING

    # Arguments:
    - `information`: Struct containting all information for model.
    - `t`: Time [s]. 
    """

    dof = information.DOF
    K = zeros(information.NX * dof, information.NX * dof)
    M = zeros(information.NX * dof, information.NX * dof)
    F = zeros(information.NX * dof)

    for ELEMENT in information.ELEMENTS
        element_global_ids = Vector{Int}(undef, length(ELEMENT.NODE_IDS) * dof)
        for (i, GLOBAL_VECTOR_ID) in enumerate(ELEMENT.GLOBAL_VECTOR_IDS)
            F[GLOBAL_VECTOR_ID] += ELEMENT.LOAD[i]
        end

        ELEMENT_K = ELEMENT.K
        ELEMENT_M = ELEMENT.M

        for i = 1:size(ELEMENT.GLOBAL_MATRIX_IDS, 1)
            GLOBAL_MATRIX_ID = ELEMENT.GLOBAL_MATRIX_IDS[i, :]
            LOCAL_MATRIX_ID = ELEMENT.LOCAL_MATRIX_IDS[i, :]

            i_L = LOCAL_MATRIX_ID[1]
            j_L = LOCAL_MATRIX_ID[2]
            K_L = ELEMENT_K[i_L, j_L]
            M_L = ELEMENT_M[i_L, j_L]

            i_G = GLOBAL_MATRIX_ID[1]
            j_G = GLOBAL_MATRIX_ID[2]
            K[i_G, j_G] = K[j_G, i_G] = K_L
            M[i_G, j_G] = M[j_G, i_G] = M_L
        end
    end

    return (K, M, F)

end
