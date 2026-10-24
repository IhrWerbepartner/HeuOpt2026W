using MHLib: MHLib 
using Random: Random
using Graphs: SimpleGraph, add_edge!

export RCkSTInstance, RCkSTSolution

"""
    Edge

Represents an edge in a graph
"""
struct Edge
    v0::Int
    v1::Int
end

"""
    EdgeCost

 - `weight`: weight of the edge in a spanning tree
 - `resource_comsumtion`: amount of resources [1..m] this edge consumes
"""
struct EdgeCost
    weight::Int
    resource_consumtion::Vector{Int}
end

"""
    RCkSTInstance

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
    edge_count::Int
    node_resources::Matrix{Int} # node [x] has resources [x, 1..n]
    graph::SimpleGraph{Int}
    edge_costs::Dict{Edge,EdgeCost}
    # jg: TODO: test if a full adjacency matrix is faster than the indirection through the hashmap.
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
    node_resources = reshape(all_values[5:resources_end_index], n, m)
    edges = permutedims(
        reshape(
            all_values[(resources_end_index + 1):(resources_end_index + edge_count * (3 + m))],
            3 + m,
            edge_count,
        ),
    )

    # jg: TODO: figure out how to rescale the resources but avoid floating point numbers (possibly (?)) for performance reasons

    edge_costs = Dict{Edge,EdgeCost}()
    sizehint!(edge_costs, size(edges, 1))

    graph = SimpleGraph(n)

    for row in axes(edges, 1)
        edge_costs[Edge(edges[row, 1], edges[row, 2])] =  EdgeCost(edges[row, 3], edges[row, 4:end])
        add_edge!(graph, edges[row, 1], edges[row, 2])
    end
    return RCkSTInstance(file_name, n, m, k, edge_count, node_resources, graph, edge_costs)
end

# jg: test reading instances is working
inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "example_instance.txt"))

mutable struct RCkSTSolution <: MHLib.SubsetVectorSolution{Int}
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
    s::RCkSTSolution, update_obj_val::Bool=true, allow_infeasible::Bool=false
) end
function MHLib.element_added_delta_eval!(
    s::RCkSTSolution, update_obj_val::Bool=true, allow_infeasible::Bool=false
) end

