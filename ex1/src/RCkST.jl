module RCkST

using MHLib:
    MHLib,
    element_added_delta_eval!,
    element_removed_delta_eval!,
    check,
    SubsetVectorSolution,
    may_be_extendible,
    to_maximize
using Graphs: SimpleGraph, add_edge!, neighbors

export RCkSTInstance,
    RCkSTSolution,
    edge_key,
    edge_id,
    set_edges!,
    selected_vertices,
    resource_slack,
    write_solution,
    read_solution,
    Edge

"""
    Edge

Represents an edge in a graph
"""
struct Edge
    v0::Int
    v1::Int
    Edge(u::Int, v::Int) = new(min(u, v), max(u, v))
    Edge((u, v)::Tuple{Int,Int}) = new(min(u, v), max(u, v))
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
- `edge_count`: total number of edges in the graph
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

    edges_hashed::Dict{Edge,Int}
    graph::SimpleGraph{Int}
    # jg: TODO: test if a full adjacency matrix is faster than the indirection through the hashmap.

    # jg: TODO: rescaling could be done by remembering the min, and range of every dimension and transforming
    # it by (x - min) / (max - min). This would map it into the inverval [0; 1].
    # We have to remember the rescaled an original resource information though to avoid floating point issues.

    # Also note that we need to take into account both vertex resources and edge consumption to determine min and max values.

    # Maybe we want a certain distribtuion for the different resources (?). So we could transform them to all be uniform or normal?
    # Could be done by an inverse CDF, like in a QQ-Plot.

    # Maybe we can compute correlation/covariance of the resource matrix. This could reduce calculation effort by first doing
    # PCA (?), Problem: the instances do not have enough dimensions for this to be relevant probably.
end

edge_key(u::Int, v::Int) = (min(u, v), max(u, v))

"""
Returns the cost of the edge spanned by `u` and `v`.
Assumes that this edge is present in the graph and errors otherwise.
"""
function edge_cost_from_vertices(inst::RCkSTInstance, u::Int, v::Int)::EdgeCost
    @assert haskey(inst.edges_hashed, Edge(u, v))
    edge_index = inst.edges_hashed[Edge(u, v)]
    return EdgeCost(inst.edges[edge_index, 3], inst.edges[edge_index, 4:end])
end

"""
    RCkSTInstance(file_name::String)

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

    if length(all_values) < 5 # some check to asssert we actually read the file
        error("invalid instance file: $file_name")
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

    resources = permutedims(reshape(all_values[5:resources_end_index], m, n))

    edges = permutedims(
        reshape(all_values[(resources_end_index + 1):end], 3 + m, edge_count)
    )

    # jg: TODO: figure out how to rescale the resources but avoid floating point numbers (possibly (?)) for performance reasons
    # we also need to take into account how much of each resource is available,
    # so in the end the ratio of resource needed / resource available should be roughly equal for all resources.

    # Create dictionary and graph
    edges_hashed = Dict{Edge,Int}()
    sizehint!(edges_hashed, edge_count)

    graph = SimpleGraph(n)

    for row in axes(edges, 1)
        u = edges[row, 1]
        v = edges[row, 2]

        1 <= u <= n && 1 <= v <= n || error("Invalid edge endpoints")

        u != v || error("Self-loop detected")

        all(x -> x >= 0, @view edges[row, 3:end]) ||
            error("Negative edge weight or consumption")

        key = Edge(edge_key(u, v))

        haskey(edges_hashed, key) && error("Duplicate edge")

        edges_hashed[key] = row

        add_edge!(graph, u, v) || error("Could not add edge")
    end

    return RCkSTInstance(
        file_name, n, m, k, edge_count, resources, edges, edges_hashed, graph
    )
end

#edge lookup helper
function edge_id(inst::RCkSTInstance, u::Int, v::Int)
    return get(inst.edges_hashed, edge_key(u, v), nothing)
end

#MISP library requires the exact names obj_val and obj_val_valid
mutable struct RCkSTSolution <: SubsetVectorSolution{Int}
    inst::RCkSTInstance
    obj_val::Int
    obj_val_valid::Bool
    x::Vector{Int}
    sel::Int
    # jg: TODO: think of some values we need for delta evaluation

    # current diff of resource "consumed" vs "accumulated". The diff is per resource.
    resource_diff::Vector{Int}

    # edges that could be added to the solution without breaking the "tree" constraint
    neighboring_edegs::Set{Int}

    # incident edges per vertex. a count of 1 specifies a leaf.
    # this means if an edge has an endpoint in a vertex with a count of 1, it is a leaf.
    incident_edges_per_vertex::Dict{Int,Int}
end

function is_leaf_edge(s::RCkSTSolution, edge::Int)
    return s.incident_edges_per_vertex[s.inst.edges[edge, 1]] == 1 ||
           s.incident_edges_per_vertex[s.inst.edges[edge, 2]] == 1
end

"""
Create a solution object from an instance.
Initializes all the delta evalutaion specific datastructures.
"""
function RCkSTSolution(inst::RCkSTInstance)
    return RCkSTSolution(
        inst,
        0,
        false,
        collect(1:inst.edge_count),
        0,
        zeros(Int, inst.m),
        Set{Int}(),
        Dict{Int,Int}(),
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
        s.sel,
        copy(s.resource_diff),
        copy(s.neighboring_edegs),
        copy(s.incident_edges_per_vertex),
    )
end

#modifying existing solution
function Base.copy!(dest::RCkSTSolution, src::RCkSTSolution)
    dest.inst === src.inst || error("Solutions belong to different instances")

    copyto!(dest.x, src.x)

    dest.sel = src.sel
    dest.obj_val = src.obj_val
    dest.obj_val_valid = src.obj_val_valid

    return dest
end

#print solution in Julia
function Base.show(io::IO, s::RCkSTSolution)
    edge_fmt = if s.sel > 20
        "$(first(s.x[1:s.sel], 10)) ... $(last(s.x[1:s.sel], 10))"
    else
        "$(s.x[1:s.sel])"
    end
    return print(
        io,
        "RCkSTSolution(sel=",
        s.sel,
        ", edges=",
        edge_fmt,
        ", weight=",
        MHLib.obj(s),
        ", slack=",
        s.resource_diff,
        ")",
    )
end

#invalidates solution
function Base.empty!(s::RCkSTSolution)
    s.sel = 0
    MHLib.invalidate!(s)
    s.resource_diff = zeros(Int, s.inst.m)
    empty!(s.neighboring_edegs)
    empty!(s.incident_edges_per_vertex)
    return s
end

#helper to set edges that are in the solution at the beginning of the Vector x
function set_edges!(s::RCkSTSolution, ids::AbstractVector{<:Integer})
    chosen_ids = Int.(ids)

    length(unique(chosen_ids)) == length(chosen_ids) || error("Duplicate edge IDs")

    all(i -> 1 <= i <= s.inst.edge_count, chosen_ids) || error("Invalid edge ID")

    chosen = Set(chosen_ids)

    s.x = vcat(chosen_ids, [i for i in 1:s.inst.edge_count if !(i in chosen)])

    s.sel = length(chosen_ids)

    MHLib.invalidate!(s)

    return s
end

#calculates the total weight of the tree
function MHLib.calc_objective(s::RCkSTSolution)
    return sum(s.inst.edges[s.x[1:s.sel], 3])
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
function selected_is_connected(s::RCkSTSolution, verts::Set{Int})
    isempty(verts) && return false

    # Build adjacency dictionary for selected edges only
    adjacency = Dict{Int,Vector{Int}}(v => Int[] for v in verts)

    for i in 1:s.sel
        e = s.x[i]

        haskey(adjacency, s.inst.edges[e, 1]) || return false
        haskey(adjacency, s.inst.edges[e, 2]) || return false
        push!(adjacency[s.inst.edges[e, 1]], s.inst.edges[e, 2])
        push!(adjacency[s.inst.edges[e, 2]], s.inst.edges[e, 1])
    end

    # Start depth-first search
    stack = [first(verts)]
    visited = Set{Int}()
    sizehint!(visited, length(verts))

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
    length(s.x) == E || error("Wrong length of solution vector x")

    0 <= s.sel <= E || error("Invalid selected-edge count")

    all(e -> 1 <= e <= E, s.x) || error("Invalid edge ID")

    #MHLib's general checks
    invoke(MHLib.check, Tuple{MHLib.SubsetVectorSolution}, s; kwargs...)

    #k-1 selected edges
    s.sel == s.inst.k - 1 ||
        error("Wrong number of selected edges, expected: $(s.inst.k - 1), got: $(s.sel)")

    #k distinct vertices
    verts = selected_vertices(s)

    length(verts) == s.inst.k || error(
        "Wrong number of selected vertices, expected: $(s.inst.k), got: $(length(verts))",
    )

    #all vertices connected (with k-1 edges & k distinct vertices -> tree)
    selected_is_connected(s, verts) || error("Selected edges are disconnected")

    #resource constraints satisfied
    slack = resource_slack(s)

    all(x -> x >= 0, slack) || error("Resource constraints violated: slack=$slack")

    return nothing
end

#write solution in solution format
function write_solution(io::IO, s::RCkSTSolution)
    instance_name = splitext(basename(s.inst.file_name))[1]

    println(io, instance_name)
    for i in 1:s.sel
        e = s.x[i]
        println(io, "$(s.inst.edges[e, 1]) $(s.inst.edges[e, 2])")
    end
end

#write solution in solution format
function write_solution(path::String, s::RCkSTSolution)
    MHLib.check(s)

    instance_name = splitext(basename(s.inst.file_name))[1]

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
    lines = [strip(line) for line in eachline(path) if !isempty(strip(line))]

    length(lines) == inst.k || error("Wrong number of solution lines")

    name = splitext(basename(inst.file_name))[1]

    lines[1] == name || error("Solution belongs to another instance")

    ids = Int[]

    for line in lines[2:end]
        tokens = split(line)

        length(tokens) == 2 || error("Invalid solution edge")

        u, v = parse.(Int, tokens)

        e = edge_id(inst, u, v)

        isnothing(e) && error("Edge ($u,$v) does not exist")

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

#function MHLib.construct!(s::RCkSTSolution, par::Int, result::Result) end
#function MHLib.local_improve!(s::RCkSTSolution, par::Int, result::Result) end
#function MHLib.shaking!(s::RCkSTSolution, par::Int, result::Result) end

"""
Returns true if the solution can be extended by another element
"""
function MHLib.may_be_extendible(s::RCkSTSolution)
    println(stderr, "may_be_extendible()")
    @assert s.sel > 0 #TODO: check if this is always the case.

    # maximum number of edges selected
    if s.sel >= s.inst.k - 1
        return false
    end

    # check for every neighboring edge it is possible to insert it without violating
    # the resource constraint.
    for e in s.neighboring_edegs
        if haskey(s.incident_edges_per_vertex, s.inst.edges[e, 1])
            neighboring_vertex = s.inst.edges[e, 2]
        else
            neighboring_vertex = s.inst.edges[e, 1]
        end
        new_resource_diff::Vector{Int} =
            s.resource_diff + s.inst.edges[e, 4:end] -
            s.inst.resources[neighboring_vertex, :]
        if all(>=(0), new_resource_diff)
            return true
        end
    end
    return false
end

"""
Delta evalution for removing and edge from the solution.

Computes the delta for the resource difference of the current solution.
This is the difference between used resources by the selected edges,
and provided resources by the incident vertices.

Updates the neighboring edges of the solution by adding neighbors that would
now not form a cycle anymore and deleting neighbors that are now not reachable by
the newly deleted vertex.

Updates the `incident_edges_per_vertex` dict for both vertices the old edge was incident to.
Deletes the vertex from the dict if no edge is incident to it anymore.

Updates the `obj_val` if `update_obj_val` is set to `true`.
"""
function MHLib.element_removed_delta_eval!(
    s::RCkSTSolution; update_obj_val::Bool=true, allow_infeasible::Bool=false
)
    removed_edge = s.x[s.sel + 1]
    incident_edges = sum(values(s.incident_edges_per_vertex))
    @assert (s.sel + 1) * 2 == incident_edges
    println(
        stderr,
        "element_removed_delta_eval!(removed: $removed_edge vs: $(s.inst.edges[removed_edge, 1:2]), allow_infeasible: $allow_infeasible)",
    )
    @assert haskey(s.incident_edges_per_vertex, s.inst.edges[removed_edge, 1]) show(s)
    @assert haskey(s.incident_edges_per_vertex, s.inst.edges[removed_edge, 2])

    if allow_infeasible
        for v in s.inst.edges[removed_edge, 1:2]
            if get(s.incident_edges_per_vertex, v, 0) == 1
                delete!(s.incident_edges_per_vertex, v)

                # add all neighboring edges that can now reached without forming a cycle.
                old_neighbors = neighbors(s.inst.graph, v)
                for n in old_neighbors
                    edge_index = s.inst.edges_hashed[Edge(v, n)]
                    if haskey(s.incident_edges_per_vertex, n)
                        push!(s.neighboring_edegs, edge_index)
                    else
                        # remove neighboring edges that are now no longer reachable.
                        delete!(s.neighboring_edegs, edge_index)
                    end
                end
            else
                s.incident_edges_per_vertex[v] -= 1
            end
        end

        s.resource_diff -= s.inst.edges[removed_edge, 4:end]

        # delta evaluate new solution
        if update_obj_val
            s.obj_val -= s.inst.edges[removed_edge, 3] # subtract edge weight
        end
        @assert s.sel * 2 == sum(values(s.incident_edges_per_vertex))
        @assert incident_edges - 2 == sum(values(s.incident_edges_per_vertex))

        # first ensure tree constraint second ensure resource constraint
        return length(s.incident_edges_per_vertex) == s.sel + 1 &&
               all(>=(0), s.resource_diff) &&
               selected_is_connected(s, Set(keys((s.incident_edges_per_vertex))))
    end

    v1 = s.inst.edges[removed_edge, 1]
    v2 = s.inst.edges[removed_edge, 2]
    vs = [v1, v2]
    removed_vertices::Vector{Int} = []

    new_resource_diff = s.resource_diff + s.inst.edges[removed_edge, 4:end]
    for v in vs
        if get(s.incident_edges_per_vertex, v, 0) == 1
            new_resource_diff -= s.inst.resources[v, :]
            delete!(s.incident_edges_per_vertex, v)
            push!(removed_vertices, v)
        else
            s.incident_edges_per_vertex[v] -= 1
        end
    end

    if !(
        length(s.incident_edges_per_vertex) == s.sel + 1 &&
        all(>=(0), new_resource_diff) &&
        selected_is_connected(s, Set(keys((s.incident_edges_per_vertex))))
    )
        # invalid solution, revert datastructures to before we did something.
        for v in vs
            if !haskey(s.incident_edges_per_vertex, v)
                new_resource_diff += s.inst.resources[v, :]
                s.incident_edges_per_vertex[v] = 1
            else
                s.incident_edges_per_vertex[v] += 1
            end
        end
        s.sel += 1
        @assert s.sel * 2 == sum(values(s.incident_edges_per_vertex))
        @assert incident_edges == sum(values(s.incident_edges_per_vertex))
        return false
    end

    s.resource_diff = new_resource_diff

    # remove all neighboring edges that can are no longer reachable.
    # TODO: do this correct even though we had a wrong solution before.
    # (consider an edge that was not connected to anything that got removed)
    # for now this is unused, but this is certainly not correct behavior.
    for v in removed_vertices
        removed_neighbors = neighbors(s.inst.graph, v)
        for n in removed_neighbors
            edge_index = s.inst.edges_hashed[Edge(v, n)]
            if haskey(s.incident_edges_per_vertex, n)
                push!(s.neighboring_edegs, edge_index)
            else
                # remove neighboring edges that are now no longer reachable.
                delete!(s.neighboring_edegs, edge_index)
            end
        end
    end

    # delta evaluate new solution
    if update_obj_val
        s.obj_val -= s.inst.edges[removed_edge, 3] # subtract edge weight
    end

    #if 98 in removed_vertices
    #    error("removed vertex 98")
    #end

    @assert s.sel * 2 == sum(values(s.incident_edges_per_vertex))
    @assert incident_edges - 2 == sum(values(s.incident_edges_per_vertex))

    return true
end

"""
Delta evalution for adding and edge to the solution.

Computes the delta for the resource difference of the current solution.
This is the difference between used resources by the selected edges,
and provided resources by the incident vertices.

Updates the neighboring edges of the solution by deleting neighbors that would
now form a cycle and adding neighbors that are now reachable by the newly added vertex.

Updates the `incident_edges_per_vertex` dict for both vertices the new edge is incident to.

Updates the `obj_val` if `update_obj_val` is set to `true`.
"""
function MHLib.element_added_delta_eval!(
    s::RCkSTSolution; update_obj_val::Bool=true, allow_infeasible::Bool=false
)
    incident_edges = sum(values(s.incident_edges_per_vertex))
    @assert (s.sel - 1) * 2 == incident_edges
    added_edge = s.x[s.sel]
    println(
        stderr,
        "element_added_delta_eval!(added: $added_edge, vs: $(s.inst.edges[added_edge, 1:2]), allow_infeasible: $allow_infeasible)",
    )

    if allow_infeasible
        for v in s.inst.edges[added_edge, 1:2]
            if !haskey(s.incident_edges_per_vertex, v)
                s.resource_diff += s.inst.resources[v, :]
                s.incident_edges_per_vertex[v] = 1

                # add all neighboring edges that can ONLY have been reached by `v`.
                added_neighbors = neighbors(s.inst.graph, v)
                for n in added_neighbors
                    edge_index = s.inst.edges_hashed[Edge(v, n)]
                    if haskey(s.incident_edges_per_vertex, n)
                        # remove neighboring edges that would now form a cycle.
                        delete!(s.neighboring_edegs, edge_index)
                    else
                        push!(s.neighboring_edegs, edge_index)
                    end
                end
            else
                s.incident_edges_per_vertex[v] += 1
            end
            @assert s.incident_edges_per_vertex[v] > 0
        end
        s.resource_diff -= s.inst.edges[added_edge, 4:end]

        # delta evaluate new solution
        if update_obj_val
            s.obj_val += s.inst.edges[added_edge, 3] # add edge weight
        end
        @assert haskey(s.incident_edges_per_vertex, s.inst.edges[added_edge, 1])
        @assert haskey(s.incident_edges_per_vertex, s.inst.edges[added_edge, 2])
        @assert s.sel * 2 == sum(values(s.incident_edges_per_vertex))

        # first ensure tree constraint second ensure resource constraint
        return length(s.incident_edges_per_vertex) == s.sel + 1 &&
               all(>=(0), s.resource_diff) &&
               selected_is_connected(s, Set(keys((s.incident_edges_per_vertex))))
    end

    v1 = s.inst.edges[added_edge, 1]
    v2 = s.inst.edges[added_edge, 2]
    vs = [v1, v2]
    added_vertices::Vector{Int} = []

    new_resource_diff::Vector{Int} = s.resource_diff - s.inst.edges[added_edge, 4:end]
    for v in vs
        if !haskey(s.incident_edges_per_vertex, v)
            new_resource_diff += s.inst.resources[v, :]
            s.incident_edges_per_vertex[v] = 1
            push!(added_vertices, v)
        else
            s.incident_edges_per_vertex[v] += 1
        end
    end

    if !(
        length(s.incident_edges_per_vertex) == s.sel + 1 &&
        all(>=(0), new_resource_diff) &&
        selected_is_connected(s, Set(keys((s.incident_edges_per_vertex))))
    )
        for v in vs
            if s.incident_edges_per_vertex[v] == 1
                delete!(s.incident_edges_per_vertex, v)
            else
                s.incident_edges_per_vertex[v] -= 1
            end
        end
        s.sel -= 1
        println(stderr, "could not add $added_edge, sel: $(s.sel)")
        @assert s.sel * 2 == sum(values(s.incident_edges_per_vertex))
        @assert incident_edges == sum(values(s.incident_edges_per_vertex))
        return false
    end

    s.resource_diff = new_resource_diff

    # add all neighboring edges that can be reached by the added vertices
    for v in added_vertices
        added_neighbors = neighbors(s.inst.graph, v)
        for n in added_neighbors
            edge_index = s.inst.edges_hashed[Edge(v, n)]
            if !(haskey(s.incident_edges_per_vertex, n))
                push!(s.neighboring_edegs, edge_index)
            else
                # remove neighboring edges that would now form a cycle.
                delete!(s.neighboring_edegs, edge_index)
            end
        end
    end

    # delta evaluate new solution
    if update_obj_val
        s.obj_val += s.inst.edges[added_edge, 3] # add edge weight
    end
    @assert haskey(s.incident_edges_per_vertex, s.inst.edges[added_edge, 1])
    @assert haskey(s.incident_edges_per_vertex, s.inst.edges[added_edge, 2])

    #if 98 in added_vertices
    #    error("added vertex 98")
    #end

    @assert s.sel * 2 == sum(values(s.incident_edges_per_vertex))
    @assert incident_edges + 2 == sum(values(s.incident_edges_per_vertex))

    return true
end

end # module