#=
This file contains the implementaion for task 1 and 2 of Exercise 1.

Task 1 description:
    Develop a meaningful deterministic construction heuristic.
Task 2 description:
    Derive a randomized construction heuristic to be applied iteratively,
    yielding diverse promising solutions.
=#

include("RCkST.jl")
using Graphs: neighbors
using Random: seed!
using MHLib: check, element_added_delta_eval!, element_removed_delta_eval!
using .RCkST:
    RCkSTInstance,
    edge_key,
    write_solution,
    Edge,
    EdgeCost,
    edge_cost_from_vertices,
    RCkSTSolution

"""
return all the edges (`v0` corresponds to a node currently in the solution, `v1` the new node)
with its corresponding cost that we could insert into the
solution given the current `selected_edges`.
"""
function neighborhood!(s::RCkSTSolution, unvisited_neighbor_index::Int)::Bool
    @assert length(s.incident_edges_per_vertex) == s.sel + 1
    original_solution_size = s.sel

    sorted_neighbors = sort(collect(s.neighboring_edegs))

    added_vertices::Vector{Int} = []
    sizehint!(added_vertices, length(sorted_neighbors))
    for n in sorted_neighbors
        if haskey(s.incident_edges_per_vertex, s.inst.edges[n, 1])
            push!(added_vertices, s.inst.edges[n, 2])
        elseif haskey(s.incident_edges_per_vertex, s.inst.edges[n, 2])
            push!(added_vertices, s.inst.edges[n, 1])
        else
            error("invalid edge to add")
        end
    end

    worst_case_resource_use = minimum(
        reshape(s.resource_diff, 1, s.inst.m) .- s.inst.edges[sorted_neighbors, 4:end] +
        s.inst.resources[added_vertices, :];
        dims=2,
    )
    edge_weights = s.inst.edges[sorted_neighbors, 3]

    # first criterion is worst case resource use, second is edge weight
    edges_sorted = sortperm(
        1:length(sorted_neighbors);
        by=x -> (worst_case_resource_use[x], -edge_weights[x]),
        rev=true,
    )

    found_neighbor::Bool = false

    for e_index in edges_sorted[unvisited_neighbor_index:end]
        if worst_case_resource_use[e_index] < 0
            break # all subsequent elements will only be worse, no need check them
        end
        s.sel += 1
        s.x[s.sel] = sorted_neighbors[e_index]
        if element_added_delta_eval!(s)
            found_neighbor = true
            @assert worst_case_resource_use[e_index] >= 0
            break
        end
        #@assert length(setdiff(sorted_neighbors, s.neighboring_edegs)) == 0
    end

    if found_neighbor
        @assert original_solution_size + 1 == s.sel
    else
        @assert original_solution_size == s.sel
    end

    return found_neighbor
end

function find_inital_vertex!(s::RCkSTSolution, sort::Vector{Int}, start_index::Int)::Int
    s.sel += 1
    s.x[s.sel] = sort[start_index]
    while !element_added_delta_eval!(s)
        @assert s.sel < s.inst.edge_count
        @assert start_index < s.inst.edge_count
        start_index += 1
        s.sel += 1
        s.x[s.sel] = sort[start_index]
    end
    return start_index
end

"""
### Greedy construction method.

    deterministc_construction_delta!(s::RCkSTSolution)::Nothing

Performs the following algorithm:
- Select an inital node based on the highest difference between resource consumption
of an edge and resources provided by its endpoints.
    - For each neighbor of all edges in the current solution calculate the above difference
    - Select the best option, for ties select the lower weighted edge.
    - If the neighbor does not fit try the next neighbors with increasing costs, meaning worse lower bounds
    - If we can not find a neighbor to add for the current solution backtrack one level
    and add the second best neighbor of the previous iteration and try again.

The idea is to create a solution that looks at the resources and tries to balance them as good
as possible, leaving much resources available for local/shaking search methods.
"""
function deterministc_construction_delta!(s::RCkSTSolution)::Nothing
    worst_case_resource_use = minimum(
        s.inst.resources[s.inst.edges[:, 1], :] +
        s.inst.resources[s.inst.edges[:, 2], :] - s.inst.edges[:, 4:end];
        dims=2,
    )
    edge_weights = s.inst.edges[:, 3]

    # first criterion is worst case resource use, second is edge weight
    initial_edge_costs = sortperm(
        1:s.inst.edge_count;
        by=x -> (worst_case_resource_use[x], -edge_weights[x]),
        rev=true,
    )

    initial_node_index = find_inital_vertex!(s, initial_edge_costs, 1)

    # we need to keep track of the best valid neighbor we currently
    # have checked if we paint ourselves into a corner.
    best_valid_neighbor::Vector{Int} = [initial_node_index]
    sizehint!(best_valid_neighbor, s.inst.k)

    while s.sel < s.inst.k - 1
        found_neighbor = false
        if neighborhood!(s, best_valid_neighbor[end])
            push!(best_valid_neighbor, 1)
            found_neighbor = true
        else
            s.sel -= 1
            element_removed_delta_eval!(s)
            if length(best_valid_neighbor) > 1
                pop!(best_valid_neighbor)
                best_valid_neighbor[end] += 1
                println(
                    stderr,
                    "backtracking, sol size: $(s.sel), diff: $(s.resource_diff), n: $best_valid_neighbor",
                )
            else
                # we need to find another initial starting positon
                initial_node_index += 1
                println(stdout, "restart, index: $initial_node_index")
                initial_node_index = find_inital_vertex!(
                    s, initial_edge_costs, initial_node_index
                )
            end
        end
        @assert s.sel == length(best_valid_neighbor)
        @debug begin
            if s.sel % 10 == 0
                println(stdout, "current solution size: $(s.sel), diff: $(s.resource_diff)")
            end
        end
    end
    s.x = [s.x[1:s.sel]; setdiff(1:s.inst.edge_count, s.x[1:s.sel])]
    s.obj_val_valid = true
    return nothing
end

module WeightedPermutation
using Random: seed!

"""
Idea taken from

- https://doi.org/10.1016/j.ipl.2005.11.003
- https://stackoverflow.com/a/15205104
"""
struct Iterator
    data
    permutation::Vector{Int}
    function Iterator(data, weights::AbstractVector{<:Real})::Iterator
        @assert length(data) <= length(weights)
        @assert all(>=(0), weights)
        w = weights[1:length(data)]
        perm = rand(length(w)) .^ (1 ./ w)
        return new(data, sortperm(perm; rev=true))
    end
    function Iterator(data, weights::AbstractVector{<:Real}, seed::Int)::Iterator
        @assert length(data) <= length(weights)
        @assert all(>=(0), weights)
        seed!(seed)
        w = weights[1:length(data)]
        perm = rand(length(w)) .^ (1 ./ w)
        return new(data, sortperm(perm; rev=true))
    end
end

function Base.iterate(iter::Iterator, state=1)
    if state > length(iter.data)
        return nothing
    end
    return (iter.data[iter.permutation[state]], state + 1)
end

"""
densitiy value of the geometrc distribtuion for `k` and probably `p`.
"""
function dgeom(k::AbstractArray, prob::Real)
    return ((1 - prob) .^ k) * prob
end
end


"""
Identical to the deterministic version, but randomize the best neighbors a bit.
This is controlled by `weights`.
"""
function randomized_neighborhood!(
    s::RCkSTSolution, unvisited_neighbor_index::Int, weights::Vector{<:Real}
)::Bool
    @assert length(s.incident_edges_per_vertex) == s.sel + 1
    original_solution_size = s.sel

    sorted_neighbors = sort(collect(s.neighboring_edegs))

    added_vertices::Vector{Int} = []
    sizehint!(added_vertices, length(sorted_neighbors))
    for n in sorted_neighbors
        if haskey(s.incident_edges_per_vertex, s.inst.edges[n, 1])
            append!(added_vertices, s.inst.edges[n, 2])
        elseif haskey(s.incident_edges_per_vertex, s.inst.edges[n, 2])
            append!(added_vertices, s.inst.edges[n, 1])
        else
            error("invalid edge to add")
        end
    end

    worst_case_resource_use = minimum(
        reshape(s.resource_diff, 1, s.inst.m) .- s.inst.edges[sorted_neighbors, 4:end] +
        s.inst.resources[added_vertices, :];
        dims=2,
    )
    edge_weights = s.inst.edges[sorted_neighbors, 3]

    # first criterion is worst case resource use, second is edge weight
    edges_sorted = sortperm(
        1:length(sorted_neighbors);
        by=x -> (worst_case_resource_use[x], -edge_weights[x]),
        rev=true,
    )

    found_neighbor::Bool = false

    random_iter = WeightedPermutation.Iterator(edges_sorted[unvisited_neighbor_index:end], weights, s.x[s.sel])

    for e_index in random_iter
        s.sel += 1
        s.x[s.sel] = sorted_neighbors[e_index]
        if element_added_delta_eval!(s)
            found_neighbor = true
            @assert worst_case_resource_use[e_index] >= 0
            break
        end
    end

    if found_neighbor
        @assert original_solution_size + 1 == s.sel
    else
        @assert original_solution_size == s.sel
    end

    return found_neighbor
end


"""
### Randomized construction method.

    randomized_construction_delta!(s::RCkSTSolution; seed=nothing)::Nothing

Idea for the randomized consturction method:

Generate a geometric distribtuion for the solution s.t. very bad solutions do not get
picked but we also do not pick the one with the least resource consumption all the time.
This could lead to picking an edge which weight is lower but resource consumption only slightly higher.
It is implemented an iterator over the array to make sure that
we can resume iteration at any point for an earlier picked edge if we land in an empty neighborhood.

The following selection of nodes is always randomized by the idea above.
- Select an inital node based on the highest difference between resource consumption
of an edge and resources provided by its endpoints.
    - For each neighbor of all edges in the current solution calculate the above difference
    - Select the best option, for ties select the lower weighted edge.
    - If the neighbor does not fit try the next neighbors with increasing costs, meaning worse lower bounds
    - If we can not find a neighbor to add for the current solution backtrack one level
    and add the second best neighbor of the previous iteration and try again.

The idea is to by chance include neighbors with better weight characteristics, then the deterministic solution.
"""
function randomized_construction_delta!(s::RCkSTSolution; seed=nothing)::Nothing
    worst_case_resource_use = minimum(
        s.inst.resources[s.inst.edges[:, 1], :] +
        s.inst.resources[s.inst.edges[:, 2], :] - s.inst.edges[:, 4:end];
        dims=2,
    )

    isnothing(seed) ? seed!() : seed!(seed)


    edge_weights = s.inst.edges[:, 3]

    # first criterion is worst case resource use, second is edge weight
    initial_edge_costs = sortperm(
        1:s.inst.edge_count;
        by=x -> (worst_case_resource_use[x], -edge_weights[x]),
        rev=true,
    )
    
    #permute the initial nodes slightly
    initial_edge_costs = initial_edge_costs[WeightedPermutation.Iterator(
        initial_edge_costs, WeightedPermutation.dgeom(1:s.inst.edge_count, 0.5)
    ).permutation]

    initial_node_index = find_inital_vertex!(s, initial_edge_costs, 1)

    random_weights = WeightedPermutation.dgeom(1:s.inst.edge_count, 0.6)

    # we need to keep track of the best valid neighbor we currently
    # have checked if we paint ourselves into a corner.
    best_valid_neighbor::Vector{Int} = [initial_node_index]
    sizehint!(best_valid_neighbor, s.inst.k)

    while s.sel < s.inst.k - 1
        found_neighbor = false
        if randomized_neighborhood!(s, best_valid_neighbor[end], random_weights)
            append!(best_valid_neighbor, 1)
            found_neighbor = true
        else
            s.sel -= 1
            element_removed_delta_eval!(s)
            if length(best_valid_neighbor) > 1
                pop!(best_valid_neighbor)
                best_valid_neighbor[end] += 1
                println(
                    stderr,
                    "backtracking, sol size: $(s.sel), diff: $(s.resource_diff), n: $(last(best_valid_neighbor, 5))",
                )
            else
                # we need to find another initial starting positon
                initial_node_index += 1
                println(stdout, "restart, index: $initial_node_index")
                initial_node_index = find_inital_vertex!(
                    s, initial_edge_costs, initial_node_index
                )
            end
        end
        @assert s.sel == length(best_valid_neighbor)
        @debug begin
            if s.sel % 10 == 0
                println(stdout, "current solution size: $(s.sel), diff: $(s.resource_diff)")
            end
        end
    end
    s.x = [s.x[1:s.sel]; setdiff(1:s.inst.edge_count, s.x[1:s.sel])]
    return nothing
end

function run_construction()
    #inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "100", "n100_k30_m2_s20019.txt"))
    #inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "200", "n200_k60_m3_s20021.txt"))
    #inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "1000", "n1000_k300_m3_s20021.txt"))
    inst = RCkSTInstance(
        joinpath(
            @__DIR__, "..", "instances", "test", "10000", "n10000_k3000_m3_s20020.txt"
        ),
    )

    rckst_sol::RCkSTSolution = RCkSTSolution(inst)
    #deterministc_construction_delta!(rckst_sol)
    randomized_construction_delta!(rckst_sol)

    check(rckst_sol)
    write_solution("sol.txt", rckst_sol)
    return println(stderr, "spanning tree weight: $(rckst_sol.obj_val)")
    #println(stderr, "slack: $(resource_difference(inst, sol))")
end

run_construction()