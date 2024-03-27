module GeneralModule

struct COORDINATES_1DC
    ID::Int
    x::Float64
end

struct COORDINATES_2DC
    ID::Int
    x::Float64
    y::Float64
end

struct COORDINATES_3DC
    ID::Int
    x::Float64
    y::Float64
    z::Float64
end

struct COORDINATES_2DB
    ID::Int
    x::Float64
    y::Float64
    theta1::Float64
end

struct MATERIAL
    ID::Int
    E::Float64
    nu::Float64
    rho::Float64
end

function fcnFromString(s)
    """
        fcnFromString(s)

    DOCSTRING

    # Arguments:
    - `s`: Function in form of string. 
    """
    f = eval(Meta.parse("t -> " * s))
    return t -> Base.invokelatest(f, t)
end

Base.@kwdef struct POINT_LOAD
    NODE_ID::Int
    DIRECTION::Int
    TIME_DEPENDENT::Bool
    MAGNITUDE::Function
end

struct DISPLACEMENT_CONDITION
    NODE_ID::Int
    DIRECTION::Int
    MAGNITUDE::Int
end

struct BOUNDARY_CONDITIONS
    POINT_LOADS::Vector{POINT_LOAD}
    DISPLACEMENT_CONDITIONS::Vector{DISPLACEMENT_CONDITION}
end

struct INFORMATION
    NX::Int
    NELX::Int
    DOF::Int
    NODAL_COORDS::Vector{Any}
    ELEMENTS::Vector{Any}
    MATERIALS::Vector{MATERIAL}
    BOUNDARY_CONDITIONS::BOUNDARY_CONDITIONS
    METHOD::String
    METHOD_ARGS::Vector{Any}
    BC_METHOD::String
end

function GET_GLOBAL_ELEMENT_IDS(NODE_IDS, DOF)
    """
        GET_GLOBAL_ELEMENT_IDS(NODE_IDS, DOF)

    DOCSTRING

    # Arguments:
    - `NODE_IDS`: Vector of node ids for a given element. 
    - `DOF`: Degrees of freedom.
    """
    N_NODES = length(NODE_IDS)
    N_IDS = N_NODES * DOF
    UPPER_TRIANGLE_N_TERMS = cumsum(1:N_IDS)[end]
    GLOBAL_VECTOR_IDS = Vector{Int}(undef, N_IDS)
    GLOBAL_MATRIX_IDS = Matrix{Int}(undef, UPPER_TRIANGLE_N_TERMS, 2)

    for (i, NODE_ID) in enumerate(NODE_IDS)
        GLOBAL_VECTOR_IDS[(i-1)*DOF+1:i*DOF] = (NODE_ID-1)*DOF+1:NODE_ID*DOF
    end
    IDX_OFFSETS = append!([0], cumsum(reverse(1:N_IDS)))
    for (i, VECTOR_GLOBAL_ID) in enumerate(GLOBAL_VECTOR_IDS)
        for (j, LOCAL_ID) in enumerate(1+i-1:N_IDS)
            MATRIX_ID = j + IDX_OFFSETS[i]
            GLOBAL_MATRIX_IDS[MATRIX_ID, :] =
                [GLOBAL_VECTOR_IDS[i] GLOBAL_VECTOR_IDS[LOCAL_ID]]
        end
    end

    return (GLOBAL_VECTOR_IDS, GLOBAL_MATRIX_IDS)

end

function GET_LOCAL_ELEMENT_IDS(NODE_IDS, DOF)
    """
        GET_LOCAL_ELEMENT_IDS(NODE_IDS, DOF)

    DOCSTRING

    # Arguments:
    - `NODE_IDS`: Vector of node ids for a given element. 
    - `DOF`: Degrees of freedom.
    """
    N_NODES = length(NODE_IDS)
    N_IDS = N_NODES * DOF
    UPPER_TRIANGLE_N_TERMS = cumsum(1:N_IDS)[end]
    LOCAL_VECTOR_IDS = 1:N_IDS
    LOCAL_MATRIX_IDS = Matrix{Int}(undef, UPPER_TRIANGLE_N_TERMS, 2)

    IDX_OFFSETS = append!([0], cumsum(reverse(1:N_IDS)))
    for (i, VECTOR_LOCAL_ID) in enumerate(LOCAL_VECTOR_IDS)
        for (j, LOCAL_ID) in enumerate(1+i-1:N_IDS)
            MATRIX_ID = j + IDX_OFFSETS[i]
            LOCAL_MATRIX_IDS[MATRIX_ID, :] =
                [LOCAL_VECTOR_IDS[i] LOCAL_VECTOR_IDS[LOCAL_ID]]
        end
    end

    return (LOCAL_VECTOR_IDS, LOCAL_MATRIX_IDS)

end

function UPDATE_1DC_COORDINATES(NODAL_COORDS, yaml)
    """
        UPDATE_1DC_COORDINATES(NODAL_COORDS, yaml)

    DOCSTRING

    # Arguments:
    - `NODAL_COORDS`: Vector to store nodal coordinates. 
    - `yaml`: YAML file in dict form using yaml package.
    """
    for NODAL_COORD in yaml["NODAL_COORDS"]
        NODAL_COORD = split(NODAL_COORD, ",")
        NODAL_COORDS[parse(Int, NODAL_COORD[1])] =
            COORDINATES_1DC(parse(Int, NODAL_COORD[1]), parse(Float64, NODAL_COORD[2]))
    end
end

function UPDATE_2DC_COORDINATES(NODAL_COORDS, yaml)
    """
        UPDATE_2DC_COORDINATES(NODAL_COORDS, yaml)

    DOCSTRING

    # Arguments:
    - `NODAL_COORDS`: Vector to store nodal coordinates. 
    - `yaml`: YAML file in dict form using yaml package.
    """
    for NODAL_COORD in yaml["NODAL_COORDS"]
        NODAL_COORD = split(NODAL_COORD, ",")
        NODAL_COORDS[parse(Int, NODAL_COORD[1])] = COORDINATES_2DC(
            parse(Int, NODAL_COORD[1]),
            parse(Float64, NODAL_COORD[2]),
            parse(Float64, NODAL_COORD[3]),
        )
    end
end

function UPDATE_3DC_COORDINATES(NODAL_COORDS, yaml)
    """
        UPDATE_3DC_COORDINATES(NODAL_COORDS, yaml)

    DOCSTRING

    # Arguments:
    - `NODAL_COORDS`: Vector to store nodal coordinates. 
    - `yaml`: YAML file in dict form using yaml package.
    """
    for NODAL_COORD in yaml["NODAL_COORDS"]
        NODAL_COORD = split(NODAL_COORD, ",")
        NODAL_COORDS[parse(Int, NODAL_COORD[1])] = COORDINATES_3DC(
            parse(Int, NODAL_COORD[1]),
            parse(Float64, NODAL_COORD[2]),
            parse(Float64, NODAL_COORD[3]),
            parse(Float64, NODAL_COORD[4]),
        )
    end
end

function UPDATE_2DB_COORDINATES(NODAL_COORDS, yaml)
    """
        UPDATE_2DB_COORDINATES(NODAL_COORDS, yaml)

    DOCSTRING

    # Arguments:
    - `NODAL_COORDS`: Vector to store nodal coordinates. 
    - `yaml`: YAML file in dict form using yaml package.
    """
    for NODAL_COORD in yaml["NODAL_COORDS"]
        NODAL_COORD = split(NODAL_COORD, ",")
        NODAL_COORDS[parse(Int, NODAL_COORD[1])] = COORDINATES_2DB(
            parse(Int, NODAL_COORD[1]),
            parse(Float64, NODAL_COORD[2]),
            parse(Float64, NODAL_COORD[3]),
        )
    end
end

end
