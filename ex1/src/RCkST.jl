using MHLib

export RCkSTInstance, RCkSTSolution

"""
- `n`: number of vertices.
- `m`: number of resources per vertex.
- `k`: edges in spanning tree.
- `resoucres`: resources available at each vertex resources[v].
- `edges`: edges in the graph with corresponding costs [1..m].
- `edges_hashed`: index of edge tuple `(v0, v1)` in `edges`. 
"""
struct RCkSTInstance
    file_name::String
    n::Int
    m::Int
    k::Int
    resources::Matrix{Int}
    edges::Matrix{Int}
    edges_hashed::Dict{Tuple{Int,Int},Int}
    # jg: TODO: test if a full adjacency matrix is faster than the indirection through the hashmap.
end

struct Edge
    v0::Int
    v1::Int
end

"""
Read RCkSTInstance from file
"""
function RCkSTInstance(file_name::String)
    all_values = Int[]
    open(file_name) do f
        for line in eachline(f)
            if !startswith(line, "#")
                for word in split(line)
                    push!(all_values, parse(Int, word))
                end
            end
        end
    end

    n = all_values[1]
    edge_count = all_values[2]
    k = all_values[3]
    m = all_values[4]

    expected_values::Int = 4 + (n * m) + (edge_count * (3 + m))
    if length(all_values) != expected_values
        error(
            "invalid number of values in RCkST instance: $(length(all_values)), expected: $expected_values",
        )
    end

    resources_end_index::Int = 5 + n * m - 1
    resources = reshape(all_values[5:resources_end_index], n, m)
    edges = permutedims(
        reshape(
            all_values[(resources_end_index + 1):(resources_end_index + edge_count * (3 + m))],
            3 + m,
            edge_count,
        ),
    )

    # jg: TODO: figure out how to rescale the resources but avoid floating point numbers (possibly (?)) for performance reasons

    index_from_edge = Dict{Tuple{Int,Int},Int}()
    sizehint!(index_from_edge, size(edges, 1))

    for row in axes(edges, 1)
        key = (edges[row, 1], edges[row, 2])
        index_from_edge[key] = row
    end
    return RCkSTInstance(file_name, n, m, k, resources, edges, index_from_edge)
end

# jg: test reading instances is working
inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "example_instance.txt"))

mutable struct RCkSTSolution <: SubsetVectorSolution{Int}
    inst::RCkSTInstance
    objective_val::Int
    objective_val_valid::Bool
    x::Vector{Int}
    sel::Int
    # jg: TODO: think of some values we need for delta evaluation

end

"""
TODO
jg: not sure what functions we actually will need, for now this is copied from the MISP example.
"""
function Base.copy(s::RCkSTSolution) end
function Base.show(io::IO, s::RCkSTSolution) end
function Base.empty!(s::RCkSTSolution) end
function MHLib.check(solution::RCkSTSolution) end
function MHLib.calc_objective(s::RCkSTSolution) end
function MHLib.construct!(s::RCkSTSolution, par::Int, result::Result) end
function MHLib.local_improve!(s::RCkSTSolution, par::Int, result::Result) end
function MHLib.shaking!(s::RCkSTSolution, par::Int, result::Result) end
function MHLib.may_be_extendible(s::RCkSTSolution) end
function MHLib.element_removed_delta_eval!(
    s::MISPSolution, update_obj_val::Bool=true, allow_infeasible::Bool=false
) end
function MHLib.element_added_delta_eval!(
    s::MISPSolution, update_obj_val::Bool=true, allow_infeasible::Bool=false
) end

