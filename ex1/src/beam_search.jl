#=
This file contains the implementaion for task 3 of Exercise 1.

Task 3 description:
    Develop a Pilot or Beam Search based heuristic.

=#

using Graphs: neighbors
using .RCkST:
    RCkSTInstance,
    Edge

struct BeamState
    selected_edges::Vector{Int}
    selected_vertices::Set{Int}
    resource_slack::Vector{Int}
    weight::Int
end

# Find all possible extensions of one partial tree
function expand_state(inst::RCkSTInstance, state::BeamState)
    successors = BeamState[]

    if isempty(state.selected_vertices)

        for e in 1:inst.edge_count

            u = inst.edges[e, 1]
            v = inst.edges[e, 2]

            new_vertices = Set([u, v])

            new_slack = [
                inst.resources[u, d] +
                inst.resources[v, d] -
                inst.edges[e, 3 + d]
                for d in 1:inst.m
            ]

            new_weight = inst.edges[e, 3]

            new_state = BeamState(
                [e],
                new_vertices,
                new_slack,
                new_weight
            )

            push!(successors, new_state)
        end

        return successors

    end
        
    for u in sort!(collect(state.selected_vertices))

        for v in neighbors(inst.graph, u)

            # Skip vertices already in the tree -> circle
            if v in state.selected_vertices
                continue
            end

            e = inst.edges_hashed[Edge(u, v)]

            new_edges = copy(state.selected_edges)
            push!(new_edges, e)

            new_vertices = copy(state.selected_vertices)
            push!(new_vertices, v)

            new_slack = copy(state.resource_slack)

            for d in 1:inst.m
                new_slack[d] += inst.resources[v, d]
                new_slack[d] -= inst.edges[e, 3 + d]
            end

            new_weight = state.weight + inst.edges[e, 3]

            push!(
                successors,
                BeamState(
                    new_edges,
                    new_vertices,
                    new_slack,
                    new_weight
                )
            )
        end
    end

    return successors
end



#prefer solutions with lower slack, if equal -> prefer lower edge weight
function f(state::BeamState)
    return (
        -minimum(state.resource_slack),
        state.weight
    )
end

# Beam Search

function beam_search(inst::RCkSTInstance, beta::Int)

    #init: x_best with nothing and S with an empty BeamState
    x_best = nothing

    S = [
        BeamState(
            Int[],
            Set{Int}(),
            zeros(Int, inst.m),
            0
        )
    ]

    while !isempty(S)

        # All candidate extensions
        C = BeamState[]

        for x in S
            #if complete solution
            if length(x.selected_vertices) == inst.k
                
                if all(d -> d >= 0, x.resource_slack)

                    if x_best === nothing ||
                       x.weight < x_best.weight

                        x_best = x
                    end
                end

            else
                # Create all one-step extensions
                append!(C, expand_state(inst, x))
            end
        end

        #remove duplicate partial trees
        seen = Set{Tuple{Vararg{Int}}}()
        distinct = BeamState[]

        for x in C
            key = Tuple(sort(x.selected_edges))

            if !(key in seen)
                push!(seen, key)
                push!(distinct, x)
            end
        end

        C = distinct

        # Keep only the best beta partial trees
        if length(C) > beta

            sort!(C, by = f)
            resize!(C, beta)

        end

        S = C
    end

    return x_best
end
