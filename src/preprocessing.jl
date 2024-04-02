using YAML
using Parameters

include("general_modules.jl")
include("element_modules.jl")
include("stiffness_modules.jl")
include("mass_module.jl")

"""
Currently code is only implmented for non-linear with 
time-dependent load.
"""

function ReadYAML(YAML_filename)
    """
        ReadYAML(YAML_filename)

    DOCSTRING

    # Arguments:
    - `YAML_filename`: YAML file name in the form of path_to_file/filename.
    """
    yaml = YAML.load(open(YAML_filename))
    NX = yaml["NX"]
    NELX = yaml["NELX"]
    METHOD = yaml["METHOD"]
    METHOD_ARGS = yaml["METHOD_ARGS"]
    if METHOD_ARGS == "nothing"
        METHOD_ARGS = [nothing]
    else
        METHOD_ARGS = parse.(Float64, split(METHOD_ARGS, ","))
    end

    BC_METHOD = yaml["BC_METHOD"]
    COORDS = yaml["COORDS"]

    if !haskey(yaml, "MASS_MATRIX")
        MASS_MATRIX = false
    else
        MASS_MATRIX = yaml["MASS_MATRIX"]
    end

    if !haskey(yaml, "DAMPING_MATRIX")
        DAMPING_MATRIX = false
    else
        DAMPING_MATRIX = yaml["DAMPING_MATRIX"]
    end

    POINT_LOADS = Vector{GeneralModule.POINT_LOAD}(
        undef,
        length(yaml["BOUNDARY_CONDITIONS"]["POINT_LOADS"]),
    )

    DISPLACEMENT_CONDITIONS = Vector{GeneralModule.DISPLACEMENT_CONDITION}(
        undef,
        length(yaml["BOUNDARY_CONDITIONS"]["DISPLACEMENT_CONDITIONS"]),
    )

    BC_DATA = yaml["BOUNDARY_CONDITIONS"]
    for (i, POINT_LOAD_DATA) in enumerate(BC_DATA["POINT_LOADS"])
        POINT_LOAD_DATA = split(POINT_LOAD_DATA, ",")
        if parse(Bool, POINT_LOAD_DATA[3])
            fcnString = POINT_LOAD_DATA[4]
            POINT_LOADS[i] = GeneralModule.POINT_LOAD(
                parse(Int, POINT_LOAD_DATA[1]),
                parse(Int, POINT_LOAD_DATA[2]),
                parse(Bool, POINT_LOAD_DATA[3]),
                GeneralModule.fcnFromString(fcnString),
            )
        else
            loadMagnitude = parse(Float64, POINT_LOAD_DATA[4])
            POINT_LOADS[i] = GeneralModule.POINT_LOAD(
                parse(Int, POINT_LOAD_DATA[1]),
                parse(Int, POINT_LOAD_DATA[2]),
                parse(Bool, POINT_LOAD_DATA[3]),
                t -> loadMagnitude,
            )
        end
    end

    for (i, DC_DATA) in enumerate(BC_DATA["DISPLACEMENT_CONDITIONS"])
        DC_DATA = split(DC_DATA, ",")
        DISPLACEMENT_CONDITIONS[i] = GeneralModule.DISPLACEMENT_CONDITION(
            parse(Int, DC_DATA[1]),
            parse(Int, DC_DATA[2]),
            parse(Float64, DC_DATA[3]),
        )
    end

    BCS = GeneralModule.BOUNDARY_CONDITIONS(POINT_LOADS, DISPLACEMENT_CONDITIONS)

    MATERIALS = Vector{GeneralModule.MATERIAL}(undef, length(yaml["MATERIALS"]))

    for (i, MATERIAL_DATA) in enumerate(yaml["MATERIALS"])
        MATERIAL_DATA = split(MATERIAL_DATA, ",")
        MATERIALS[i] = GeneralModule.MATERIAL(
            parse(Int, MATERIAL_DATA[1]),
            parse(Float64, MATERIAL_DATA[2]),
            parse(Float64, MATERIAL_DATA[3]),
            parse(Float64, MATERIAL_DATA[4]),
        )
    end

    #Currently only works with all nodes in same coordinate assumption. Can
    #change at later date!

    if COORDS == "1DC"
        NODAL_COORDS = Vector{GeneralModule.COORDINATES_1DC}(undef, NX)
        GeneralModule.UPDATE_1DC_COORDINATES(NODAL_COORDS, yaml)
        DOF = 1
    elseif COORDS == "2DC"
        GeneralModule.NODAL_COORDS = Vector{GeneralModule.COORDINATES_2DC}(undef, NX)
        UPDATE_2DC_COORDINATES(NODAL_COORDS, yaml)
        DOF = 2
    elseif COORDS == "3DC"
        NODAL_COORDS = Vector{GeneralModule.COORDINATES_3DC}(undef, NX)
        GeneralModule.UPDATE_3DC_COORDINATES(NODAL_COORDS, yaml)
        DOF = 3
    elseif COORDS == "2DB"
        NODAL_COORDS = Vector{GeneralModule.COORDINATES_1DB}(undef, NX)
        GeneralModule.UPDATE_2DB_COORDINATES(NODAL_COORDS, yaml)
        DOF = 3
    end


    ELEMENTS = Vector{Any}(undef, NELX)
    for (key, ELEMENTS_DATA) in yaml["ELEMENTS"]
        if key == "T2D2"
            for ELEMENT_DATA in ELEMENTS_DATA
                ELEMENT_DATA = split(ELEMENT_DATA, ",")
                ID = parse(Int, ELEMENT_DATA[1])
                NODE_IDS = [parse(Int, i) for i in ELEMENT_DATA[2:3]]
                LOCAL_VECTOR_IDS, LOCAL_MATRIX_IDS =
                    GeneralModule.GET_LOCAL_ELEMENT_IDS(NODE_IDS, DOF)
                GLOBAL_VECTOR_IDS, GLOBAL_MATRIX_IDS =
                    GeneralModule.GET_GLOBAL_ELEMENT_IDS(NODE_IDS, DOF)
                AREA = parse(Float64, ELEMENT_DATA[4])
                MATERIAL_ID = parse(Int, ELEMENT_DATA[5])
                Le, k, f = StiffnessModule.T2D2_STIFFNESS(
                    MATERIALS[MATERIAL_ID],
                    AREA,
                    NODAL_COORDS[NODE_IDS[1]],
                    NODAL_COORDS[NODE_IDS[2]],
                )
                if MASS_MATRIX == "Consistent"
                    m = MassModule.T2D2_CONSISTENT_MASS(
                        MATERIALS[MATERIAL_ID],
                        AREA,
                        NODAL_COORDS[NODE_IDS[1]],
                        NODAL_COORDS[NODE_IDS[2]],
                    )
                elseif MASS_MATRIX == "HRZ"
                    m = MassModule.T2D2_HRZ_LUMPED_MASS(
                        MATERIALS[MATERIAL_ID],
                        AREA,
                        NODAL_COORDS[NODE_IDS[1]],
                        NODAL_COORDS[NODE_IDS[2]],
                    )
                elseif MASS_MATRIX == "Optimal"
                    m = MassModule.T2D2_OPTIMAL_LUMPED_MASS(
                        MATERIALS[MATERIAL_ID],
                        AREA,
                        NODAL_COORDS[NODE_IDS[1]],
                        NODAL_COORDS[NODE_IDS[2]],
                    )
                else
                    m = k * 0
                end
                ELEMENTS[parse(Int, ELEMENT_DATA[1])] = ElementsModule.T2D2(
                    ID,
                    NODE_IDS,
                    AREA,
                    MATERIAL_ID,
                    Le,
                    k,
                    m,
                    f,
                    LOCAL_VECTOR_IDS,
                    LOCAL_MATRIX_IDS,
                    GLOBAL_VECTOR_IDS,
                    GLOBAL_MATRIX_IDS,
                )
            end

        elseif key == "B11"
            for ELEMENT_DATA in ELEMENTS_DATA
                ELEMENT_DATA = split(ELEMENT_DATA, ",")
                ID = parse(Int, ELEMENT_DATA[1])
                NODE_IDS = [parse(Int, i) for i in ELEMENT_DATA[2:3]]
                LOCAL_VECTOR_IDS, LOCAL_MATRIX_IDS =
                    GeneralModule.GET_LOCAL_ELEMENT_IDS(NODE_IDS, DOF)
                GLOBAL_VECTOR_IDS, GLOBAL_MATRIX_IDS =
                    GeneralModule.GET_GLOBAL_ELEMENT_IDS(NODE_IDS, DOF)
                AREA = parse(Float64, ELEMENT_DATA[4])
                MATERIAL_ID = parse(Int, ELEMENT_DATA[5])
                Le, k, f = StiffnessModule.B11_STIFFNESS(
                    MATERIALS[MATERIAL_ID],
                    AREA,
                    NODAL_COORDS[NODE_IDS[1]],
                    NODAL_COORDS[NODE_IDS[2]],
                )
                if MASS_MATRIX == "Consistent"
                    m = MassModule.B11_CONSISTENT_MASS(
                        MATERIALS[MATERIAL_ID],
                        AREA,
                        NODAL_COORDS[NODE_IDS[1]],
                        NODAL_COORDS[NODE_IDS[2]],
                    )
                elseif MASS_MATRIX == "HRZ"
                    m = MassModule.B11_HRZ_LUMPED_MASS(
                        MATERIALS[MATERIAL_ID],
                        AREA,
                        NODAL_COORDS[NODE_IDS[1]],
                        NODAL_COORDS[NODE_IDS[2]],
                    )
                elseif MASS_MATRIX == "Optimal"
                    m = MassModule.B11_OPTIMAL_LUMPED_MASS(
                        MATERIALS[MATERIAL_ID],
                        AREA,
                        NODAL_COORDS[NODE_IDS[1]],
                        NODAL_COORDS[NODE_IDS[2]],
                    )
                else
                    m = k * 0
                end

                ELEMENTS[parse(Int, ELEMENT_DATA[1])] = ElementsModule.B11(
                    ID,
                    NODE_IDS,
                    AREA,
                    MATERIAL_ID,
                    Le,
                    k,
                    m,
                    f,
                    LOCAL_VECTOR_IDS,
                    LOCAL_MATRIX_IDS,
                    GLOBAL_VECTOR_IDS,
                    GLOBAL_MATRIX_IDS,
                )
            end
        end
    end

    information = GeneralModule.INFORMATION(
        NX,
        NELX,
        DOF,
        NODAL_COORDS,
        ELEMENTS,
        MATERIALS,
        BCS,
        METHOD,
        METHOD_ARGS,
        BC_METHOD,
    )

    return information
end
