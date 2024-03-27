# GENERAL INFORMATION RELATING TO NAMING CONVENTION:
#
# Tetrahedral: TaDb, a=dimension, b=number of nodes
# Beam: Bab, a=dimension, b=element order

module ElementsModule

struct ELEMENT
    ID::Int
    NODE_IDS::Vector{Int}
    AREA::Float64
    MATERIAL_ID::Int
    LE::Float64
    K::Matrix{Float64}
    M::Matrix{Float64}
    LOAD::Vector{Float64}
    LOCAL_VECTOR_IDS::Vector{Int}
    LOCAL_MATRIX_IDS::Array{Int}
    GLOBAL_VECTOR_IDS::Vector{Int}
    GLOBAL_MATRIX_IDS::Array{Int}
end

struct T2D2
    ID::Int
    NODE_IDS::Vector{Int}
    AREA::Float64
    MATERIAL_ID::Int
    LE::Float64
    K::Matrix{Float64}
    M::Matrix{Float64}
    LOAD::Vector{Float64}
    LOCAL_VECTOR_IDS::Vector{Int}
    LOCAL_MATRIX_IDS::Array{Int}
    GLOBAL_VECTOR_IDS::Vector{Int}
    GLOBAL_MATRIX_IDS::Array{Int}
end

struct B11
    ID::Int
    NODE_IDS::Vector{Int}
    AREA::Float64
    MATERIAL_ID::Int
    LE::Float64
    K::Matrix{Float64}
    M::Matrix{Float64}
    LOAD::Vector{Float64}
    LOCAL_VECTOR_IDS::Vector{Int}
    LOCAL_MATRIX_IDS::Array{Int}
    GLOBAL_VECTOR_IDS::Vector{Int}
    GLOBAL_MATRIX_IDS::Array{Int}
end

end
