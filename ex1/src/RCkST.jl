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
- `k`: vertices in spanning tree.
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

    edge_costs = Dict{Edge,EdgeCost}()
    sizehint!(edge_costs, size(edges, 1))

    graph = SimpleGraph(n)

    # since the graph is undirected u,v == v,u -> so that they identify to the same edge
    for row in axes(edges, 1)
        key = (edges[row, 1], edges[row, 2])
        index_from_edge[key] = row
    end
    return RCkSTInstance(file_name, n, m, k, edge_count, node_resources, graph, edge_costs)
end

# jg: test reading instances is working
inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "example_instance.txt"))

#edge lookup
function edge_id(inst::RCkSTInstance, u::Int, v::Int)
    return get(
        inst.edges_hashed,
        edge_key(u, v),
        nothing
    )
end

#MISP library requires the exact names obj_val and obj_val_invalid
mutable struct RCkSTSolution <: MHLib.SubsetVectorSolution{Int}
    inst::RCkSTInstance
    obj_val::Int
    obj_val_valid::Bool
    x::Vector{Int}
    sel::Int
    # jg: TODO: think of some values we need for delta evaluation

end

function RCkSTSolution(inst::RCkSTInstance)
    return RCkSTSolution(
        inst,
        0,
        false,
        collect(1:size(inst.edges, 1)),
        0
    )
end
# RC-k-ST is a minimization problem
MHLib.to_maximize(::Type{RCkSTSolution}) = false

function Base.copy(s::RCkSTSolution)
    return RCkSTSolution(
        s.inst,
        s.obj_val,
        s.obj_val_valid,
        copy(s.x),
        s.sel
    )
end

function Base.copy!(dest::RCkSTSolution, src::RCkSTSolution)

    dest.inst === src.inst ||
        error("Solutions belong to different instances")

    copyto!(dest.x, src.x)

    dest.sel = src.sel
    dest.obj_val = src.obj_val
    dest.obj_val_valid = src.obj_val_valid

    return dest
end

function Base.show(io::IO, s::RCkSTSolution)
    print(
        io,
        "RCkSTSolution(edges=",
        s.x[1:s.sel],
        ", weight=",
        MHLib.obj(s),
        ")"
    )
end

function Base.empty!(s::RCkSTSolution)
    s.sel = 0
    MHLib.invalidate!(s)
    return s
end

#calculates the total weight of the tree
function MHLib.calc_objective(s::RCkSTSolution)
    total = 0

    for i in 1:s.sel
        e = s.x[i]
        total += s.inst.edges[e, 3]
    end

    return total
end

#helper: find selected vertices in solution
function selected_vertices(s::RCkSTSolution)

    vertices = Set{Int}()

    for i in 1:s.sel

        e = s.x[i]

        u = s.inst.edges[e, 1]
        v = s.inst.edges[e, 2]

        push!(vertices, u)
        push!(vertices, v)
    end

    return vertices
end

# Calculates the (resources provided - resources= consumed for a solution
function resource_slack(s::RCkSTSolution)

    slack = zeros(Int, s.inst.m)

    #resources provided by selected vertices
    for v in selected_vertices(s)
        for d in 1:s.inst.m
            slack[d] += s.inst.resources[v, d]
        end
    end

    #resources consumed by selected edges
    for i in 1:s.sel

        e = s.x[i]

        for d in 1:s.inst.m
            slack[d] -= s.inst.edges[e, 3 + d]
        end
    end

    return slack
end

#reading solution from a file
function read_solution(inst::RCkSTInstance, path::String)

    lines = [
        strip(line)
        for line in eachline(path)
        if !isempty(strip(line))
    ]

    length(lines) == inst.k ||
        error("Wrong number of solution lines")

    name = splitext(basename(inst.file_name))[1]

    lines[1] == name ||
        error("Solution belongs to another instance")

    ids = Int[]

    for line in lines[2:end]

        tokens = split(line)

        length(tokens) == 2 ||
            error("Invalid solution edge")

        u, v = parse.(Int, tokens)

        e = edge_id(inst, u, v)

        isnothing(e) &&
            error("Edge ($u,$v) does not exist")

        push!(ids, e)
    end

    s = set_edges!(RCkSTSolution(inst), ids)

    MHLib.check(s)

    return s
end



"""
TODO
jg: not sure what functions we actually will need, for now this is copied from the MISP example.
"""
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

