#=
This file contains the implementaion for task 1 and 2 of Exercise 1.

Task 1 description:
    Develop a meaningful deterministic construction heuristic.
Task 2 description:
    Derive a randomized construction heuristic to be applied iteratively,
    yielding diverse promising solutions.
=#

include("RCkST.jl")
using Random: seed!
using Graphs: neighbors
using MHLib: check
using .RCkST:
    RCkSTInstance,
    edge_key,
    write_solution,
    Edge,
    EdgeCost,
    edge_cost_from_vertices,
    RCkSTSolution

struct EdgeRating
    resource_diff::Int
    edge_weight::Int
end

function Base.isless(a::EdgeRating, b::EdgeRating)::Bool
    if a.resource_diff < b.resource_diff
        return true
    end
    if a.resource_diff > b.resource_diff
        return false
    end
    return a.edge_weight < b.edge_weight
end

"""
Rates a neighbor by the difference of resources the new vertex aquires
minus the cost of the edge to said vertex.
"""
function rate_neigbor(inst::RCkSTInstance, src::Int, neighbor::Int)::EdgeRating
    edge_cost::EdgeCost = edge_cost_from_vertices(inst, src, neighbor)
    return EdgeRating(
        sum(edge_cost.resource_consumtion) - sum(inst.resources[neighbor]),
        edge_cost.weight,
    )
end

"""
return all the edges (`v0` corresponds to a node currently in the solution, `v1` the new node)
with its corresponding cost that we could insert into the
solution given the current `selected_edges`.
"""
function neighborhood(
    inst::RCkSTInstance, selected_edges::Vector{Int}
)::Dict{Tuple{Int,Int},EdgeRating}
    res::Dict{Tuple{Int,Int},EdgeRating} = Dict()

    selected_nodes::Set{Int} = Set(inst.edges[selected_edges, 1:2])
    @assert length(selected_nodes) == length(selected_edges) + 1
    for selected_node in selected_nodes
        for neighbor in neighbors(inst.graph, selected_node)

            # this means the neighbor is already in our solution
            if neighbor in selected_nodes
                continue
            end

            res[(selected_node, neighbor)] = rate_neigbor(inst, selected_node, neighbor)
        end
    end
    return res
end

"""
Returns the difference for each resource of provided resources vs consumed ones.
"""
function resource_difference(
    inst::RCkSTInstance, selected_edges::Vector{Int}
)::Vector{Int}
    selected_nodes = unique(inst.edges[selected_edges, 1:2])

    diff = sum(inst.resources[selected_nodes, :]; dims=1)
    return vec(diff -= sum(inst.edges[selected_edges, 4:end]; dims=1))
end

"""
### Greedy construction method.

    deterministc_construction(inst::RCkSTInstance)::Vector{Int}

Performs the following algorithm:
- Select an inital node based on the lowest difference between resource consumption
of an edge and resources provided by its endpoints.
    - For each neighbor of all edges in the current solution calculate the above difference
    - Select the best option, for ties select the lower weighted edge.
    - If the neighbor does not fit try the next neighbors with increasing costs
    - If we can not find a neighbor to add for the current solution backtrack one level
    and add the second best neighbor of the previous iteration and try again.

Returns the indices of edges selected spannig the tree.
"""
function deterministc_construction(inst::RCkSTInstance)::Vector{Int}
    sol::Vector{Int} = []
    sizehint!(sol, inst.k - 1)
    node_resources = sum.(eachrow(inst.resources))

    initial_edge_costs = sortperm(
        sum(inst.edges[:, 4:end]) .- node_resources[inst.edges[:, 1]] -
        node_resources[inst.edges[:, 2]],
    )

    initial_node_index = 1

    append!(sol, initial_edge_costs[initial_node_index])

    # we need to keep track of the best valid neighbor we currently
    # have checked if we paint ourselves into a corner.
    best_valid_neighbor::Vector{Int} = [1]
    sizehint!(best_valid_neighbor, inst.k)

    while length(sol) < inst.k - 1
        neighbor_ratings::Vector{Pair{Tuple{Int,Int},EdgeRating}} = sort(
            collect(neighborhood(inst, sol)); by=last
        ) # node ids with corresponding ratings
        # println(stdout, "neighborhood: $neighbor_ratings")
        found_valid_neighbor = false
        for (new_edge, _) in neighbor_ratings[last(best_valid_neighbor):end]
            new_edge_index = inst.edges_hashed[Edge(new_edge)]
            @assert !(new_edge_index in sol)
            append!(sol, new_edge_index)
            if all(>=(0), resource_difference(inst, sol))
                found_valid_neighbor = true
                append!(best_valid_neighbor, 1)
                break
            else
                pop!(sol)
            end
        end
        if !found_valid_neighbor
            best_valid_neighbor[end] += 1
            println(
                "no valid neighbor found, current solution (size: $(length(sol))): $sol, best_valid_neighbor: $best_valid_neighbor",
            ) # for now
        end

        if best_valid_neighbor[end] > length(neighbor_ratings)
            if length(best_valid_neighbor) > 1
                pop!(sol)
                pop!(best_valid_neighbor)
                best_valid_neighbor[end] += 1
            else
                # we need to find another initial starting positon
                initial_node_index += 1
                sol[1] = initial_edge_costs[initial_node_index]
            end
        end
        @assert length(sol) == length(best_valid_neighbor)
        if length(sol) % 10 == 0
            println(stdout, "current solution size: $(length(sol))")
        end
    end
    return sol
end

"""
Idea taken from

- https://doi.org/10.1016/j.ipl.2005.11.003
- https://stackoverflow.com/a/15205104
"""
struct WeightedPermutation
    data
    permutation::Vector{Int}
    function WeightedPermutation(
        data, weights::AbstractVector{<:Real}
    )::WeightedPermutation
        @assert length(data) == length(weights)
        @assert all(>=(0), weights)
        perm = rand(length(weights)) .^ (1 ./ weights)
        return new(data, sortperm(perm; rev=true))
    end
    function WeightedPermutation(
        data, weights::AbstractVector{<:Real}, seed::Int
    )::WeightedPermutation
        @assert length(data) == length(weights)
        @assert all(>=(0), weights)
        seed!(seed)
        perm = rand(length(weights)) .^ (1 ./ weights)
        return new(data, sortperm(perm; rev=true))
    end
end

function Base.iterate(iter::WeightedPermutation, state=1)
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

"""
### Randomized construction method.

    randomized_construction(inst::RCkSTInstance)::Vector{Int}

Idea for the randomized consturction method:

Generate a geometric distribtuion for the solution s.t. very bad solutions do not get
picked but we also do not pick the one with the least resource consumption all the time.
This could lead to picking an edge which weight is lower but resource consumption only slightly higher.
It is implemented an iterator over the array to make sure that
we can resume iteration at any point for an earlier picked edge if we land in an empty neighborhood.

The following selection of nodes is always randomized by the idea above.
- Select an inital node based on the lowest difference between resource consumption
of an edge and resources provided by its endpoints.
    - For each neighbor of all edges in the current solution calculate the above difference
    - Select the best option, for ties select the lower weighted edge.
    - If the neighbor does not fit try the next neighbors with increasing costs
    - If we can not find a neighbor to add for the current solution backtrack one level
    and add the second best neighbor of the previous iteration and try again.

Returns the indices of edges selected spannig the tree.
"""
function randomized_construction(inst::RCkSTInstance)::Vector{Int}
    sol::Vector{Int} = []
    sizehint!(sol, inst.k - 1)
    node_resources = sum.(eachrow(inst.resources))

    initial_edge_costs = sortperm(
        sum(inst.edges[:, 4:end]) .- node_resources[inst.edges[:, 1]] -
        node_resources[inst.edges[:, 2]],
    )
    random_weights = dgeom(1:inst.edge_count, 0.7)
    initial_edge_iter = WeightedPermutation(
        initial_edge_costs, random_weights[1:length(initial_edge_costs)]
    )

    (initial_edge_index, initial_edge_state) = iterate(initial_edge_iter)

    append!(sol, initial_edge_index)

    # we need to keep track of the best valid neighbor we currently
    # have checked if we paint ourselves into a corner.
    best_valid_neighbor::Vector{Int} = [1]
    sizehint!(best_valid_neighbor, inst.k)

    while length(sol) < inst.k - 1
        neighbor_ratings::Vector{Pair{Tuple{Int,Int},EdgeRating}} = sort(
            collect(neighborhood(inst, sol)); by=last
        ) # node ids with corresponding ratings
        # println(stdout, "neighborhood: $neighbor_ratings")
        found_valid_neighbor = false
        unexplored_neighbors = neighbor_ratings[last(best_valid_neighbor):end]
        for (new_edge, _) in WeightedPermutation(
            unexplored_neighbors, random_weights[1:length(unexplored_neighbors)], sol[end]
        )
            new_edge_index = inst.edges_hashed[Edge(new_edge)]
            @assert !(new_edge_index in sol)
            append!(sol, new_edge_index)
            if all(>=(0), resource_difference(inst, sol))
                found_valid_neighbor = true
                append!(best_valid_neighbor, 1)
                break
            else
                pop!(sol)
            end
        end
        if !found_valid_neighbor
            best_valid_neighbor[end] += 1
            # println( "no valid neighbor found, current solution (size: $(length(sol))): $sol, best_valid_neighbor: $best_valid_neighbor",) # for now
        end

        if best_valid_neighbor[end] > 30 # could also be, length(neighbor_ratings) but this can lead to very long searches in resource intensive neighborhoods where we do not find a quick way out.
            if length(best_valid_neighbor) > 1
                pop!(sol)
                pop!(best_valid_neighbor)
                best_valid_neighbor[end] += 1
            else
                # we need to find another initial starting positon
                iter_result = iterate(initial_edge_iter, initial_edge_state)
                if iter_result === nothing
                    error("exhausted all starting nodes somehow") # for now
                end
                sol[1] = iter_result[1]
                initial_edge_state = iter_result[2]
            end
        end
        @assert length(sol) == length(best_valid_neighbor)
        if length(sol) % 10 == 0
            println(stdout, "current solution size: $(length(sol))")
        end
    end
    return sol
end

#inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "100", "n100_k30_m2_s20019.txt"))
#inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "200", "n200_k60_m3_s20021.txt"))
#inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "1000", "n1000_k300_m3_s20021.txt"))
inst = RCkSTInstance(
    joinpath(@__DIR__, "..", "instances", "test", "10000", "n10000_k3000_m3_s20020.txt")
)
# sol = deterministc_construction(inst)
sol = randomized_construction(inst)

# =============== PRINTING ================
rckst_sol::RCkSTSolution = RCkSTSolution(
    inst,
    sum(inst.edges[sol, 3]),
    true,
    [sol; setdiff(1:inst.edge_count, sol)],
    length(sol),
)
check(rckst_sol)
write_solution(stdout, rckst_sol)
println(stderr, "spanning tree weight: $(rckst_sol.obj_val)")