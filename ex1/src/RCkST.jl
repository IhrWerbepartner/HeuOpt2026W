import MHLib
using Graphs: SimpleGraph, add_edge!

export RCkSTInstance, RCkSTSolution,
       edge_key, edge_id, set_edges!,
       selected_vertices, resource_slack,
       write_solution, read_solution

"""
    Edge

Represents an edge in a graph

struct Edge
    v0::Int
    v1::Int
end
-> not needed
"""

"""
    EdgeCost

 - `weight`: weight of the edge in a spanning tree
 - `resource_comsumtion`: amount of resources [1..m] this edge consumes

struct EdgeCost
    weight::Int
    resource_consumtion::Vector{Int}
end
-> not needed
"""

"""
    RCkSTInstance

- `n`: number of vertices.
- `m`: number of resources per vertex.
- `k`: vertices in spanning tree.
-  `edge_count`: total number of edges in the graph
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

    resources::Matrix{Int}
    edges::Matrix{Int}

    edges_hashed::Dict{Tuple{Int,Int},Int}
    graph::SimpleGraph{Int}
    # jg: TODO: test if a full adjacency matrix is faster than the indirection through the hashmap.
end

edge_key(u::Int, v::Int) =
    (min(u, v), max(u, v))

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

    resources_end_index = 4 + n * m

    resources = permutedims(
    reshape(
        all_values[5:resources_end_index],
        m, n
    )
    )

    edges = permutedims(
        reshape(
            all_values[(resources_end_index + 1):end],
            3 + m,
            edge_count,
        ),
    )


    # jg: TODO: figure out how to rescale the resources but avoid floating point numbers (possibly (?)) for performance reasons

     # Create dictionary and graph
    edges_hashed = Dict{Tuple{Int,Int},Int}()
    sizehint!(edges_hashed, edge_count)

    graph = SimpleGraph(n)

    for row in axes(edges, 1)

        u = edges[row, 1]
        v = edges[row, 2]

        1 <= u <= n && 1 <= v <= n ||
            error("Invalid edge endpoints")

        u != v || error("Self-loop detected")

        all(x -> x >= 0, @view edges[row, 3:end]) ||
            error("Negative edge weight or consumption")

        key = edge_key(u, v)

        haskey(edges_hashed, key) &&
            error("Duplicate edge")

        edges_hashed[key] = row

        add_edge!(graph, u, v) ||
            error("Could not add edge")
    end

    return RCkSTInstance(
        file_name,
        n,
        m,
        k,
        edge_count,
        resources,
        edges,
        edges_hashed,
        graph
    )
end

# jg: test reading instances is working
inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "example_instance.txt"))

#edge lookup helper
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
        collect(1:inst.edge_count),
        0
    )
end

# RC-k-ST is a minimization problem
MHLib.to_maximize(::Type{RCkSTSolution}) = false

#creating new solutions (when previous solutions shouldnt be changed)
function Base.copy(s::RCkSTSolution)
    return RCkSTSolution(
        s.inst,
        s.obj_val,
        s.obj_val_valid,
        copy(s.x),
        s.sel
    )
end

#modifying existing solution
function Base.copy!(dest::RCkSTSolution, src::RCkSTSolution)

    dest.inst === src.inst ||
        error("Solutions belong to different instances")

    copyto!(dest.x, src.x)

    dest.sel = src.sel
    dest.obj_val = src.obj_val
    dest.obj_val_valid = src.obj_val_valid

    return dest
end

#print solution in Julia
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

#invalidates solution
function Base.empty!(s::RCkSTSolution)
    s.sel = 0
    MHLib.invalidate!(s)
    return s
end

#helper to set edges that are in the solution at the beginning of the Vector x
function set_edges!(
    s::RCkSTSolution,
    ids::AbstractVector{<:Integer}
)

    chosen_ids = Int.(ids)

    length(unique(chosen_ids)) == length(chosen_ids) ||
        error("Duplicate edge IDs")

    all(i -> 1 <= i <= s.inst.edge_count, chosen_ids) ||
        error("Invalid edge ID")

    chosen = Set(chosen_ids)

    s.x = vcat(
        chosen_ids,
        [
            i for i in 1:s.inst.edge_count
            if !(i in chosen)
        ]
    )

    s.sel = length(chosen_ids)

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

        push!(vertices, s.inst.edges[e, 1])
        push!(vertices, s.inst.edges[e, 2])
    end

    return vertices
end

# Calculates the (resources provided - resources consumed) for a solution
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

#helper: checking if selections form a tree, using DFS search
function selected_is_connected(
    s::RCkSTSolution,
    verts::Set{Int}
)

    isempty(verts) && return false

    # Build adjacency dictionary for selected edges only
    adjacency = Dict{Int,Vector{Int}}(
        v => Int[] for v in verts
    )

    for i in 1:s.sel
        e = s.x[i]

        push!(adjacency[s.inst.edges[e, 1]], s.inst.edges[e, 2])
        push!(adjacency[s.inst.edges[e, 2]], s.inst.edges[e, 1])
    end

    # Start depth-first search
    stack = [first(verts)]
    visited = Set{Int}()

    while !isempty(stack)

        v = pop!(stack)

        v in visited && continue

        push!(visited, v)

        for neighbor in adjacency[v]
            if !(neighbor in visited)
                push!(stack, neighbor)
            end
        end
    end

    return length(visited) == length(verts)
end

#checker of valid solution
function MHLib.check(s::RCkSTSolution; kwargs...)

    E = s.inst.edge_count

    #solution representation
    length(s.x) == E ||
        error("Wrong length of solution vector x")

    0 <= s.sel <= E ||
        error("Invalid selected-edge count")

    all(e -> 1 <= e <= E, s.x) ||
    error("Invalid edge ID")

    #MHLib's general checks
    invoke(
        MHLib.check,
        Tuple{MHLib.SubsetVectorSolution},
        s;
        kwargs...
    )

    #k-1 selected edges
    s.sel == s.inst.k - 1 ||
        error("Wrong number of selected edges")

    #k distinct vertices
    verts = selected_vertices(s)

    length(verts) == s.inst.k ||
        error("Wrong number of selected vertices")

    #all vertices connected (with k-1 edges & k distinct vertices -> tree)
    selected_is_connected(s, verts) ||
        error("Selected edges are disconnected")

    #resource constraints satisfied
    slack = resource_slack(s)

    all(x -> x >= 0, slack) ||
        error("Resource constraints violated: slack=$slack")

    return nothing
end

#write solution in solution format
function write_solution(path::String, s::RCkSTSolution)

    MHLib.check(s)

    instance_name =
        splitext(basename(s.inst.file_name))[1]

    open(path, "w") do io

        println(io, instance_name)

        for i in 1:s.sel

            e = s.x[i]

            println(io, "$(s.inst.edges[e, 1]) $(s.inst.edges[e, 2])")
        end
    end

    return path
end



#reading solution from a file and creating a RCkSTSolution
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
"""
