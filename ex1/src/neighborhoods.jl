include("deterministic_construction.jl")
using Graphs: neighbors
using MHLib:
    Result,
    MHMethod,
    construct!,
    local_improve!,
    shaking!,
    two_exchange_random_fill_neighborhood_search!,
    remove_randomly_selected!,
    fillup!,
    GVNS,
    run!,
    method_statistics,
    MHLib,
    main_results
using .RCkST:
    RCkSTInstance,
    edge_key,
    write_solution,
    Edge,
    EdgeCost,
    edge_cost_from_vertices,
    RCkSTSolution

"""
Local improvement can follow two paths:

1. Minimize the objective value (edge weight) while consuming all avaialbe resource budget.
2. Minimize the objective value while not consuming more resoruces than the current solution.
Furthermore this can be split by:
    - strictly not consuming more of any resource than currently
    - the sum of resource consuption (scaled) does not exceed the current solution.

"1." has the potential of finding a local optimum, but restricts further move by consuming all possible budget.

"2." may only find little improvements in terms of objective value but allows much more flexibility
in the solution which can be utilized by further moves.


A strategy for 2. could be:
- Sort leaf edges (can be removed without breaking the tree property of the solution) by cost difference between
edge consumption and leaf vertex resource gain.
- Swap the most expensive edges for neighbors ((present_vertex) --- edge --> (neigbor_vertex)) which are of lower cost
than the edge we would remove.
"""

function MHLib.construct!(s::RCkSTSolution, ::Nothing, ::Result)
    println(stderr, "here")
    deterministc_construction_delta!(s)
    println(stderr, s)
end

function MHLib.local_improve!(s::RCkSTSolution, ::Nothing, result::Result)
    println(stderr, "there")
    println(stderr, s)
    if !two_exchange_random_fill_neighborhood_search!(s, false)
        println(stderr, "no change")
        result.changed = false
    end
end

function MHLib.shaking!(s::RCkSTSolution, par::Int, ::Result)
    remove_randomly_selected!(s, par)
    return fillup!(s)
end

function run_local_search()
    inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "100", "n100_k30_m2_s20019.txt"))
    #inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "200", "n200_k60_m3_s20021.txt"))
    #inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "1000", "n1000_k300_m3_s20021.txt"))
    #inst = RCkSTInstance(
    #    joinpath(
    #        @__DIR__, "..", "instances", "test", "10000", "n10000_k3000_m3_s20020.txt"
    #    ),
    #)

    sol::RCkSTSolution = RCkSTSolution(inst)

    # we apply a variable neighborhod search
    alg = GVNS(
        sol,
        [MHMethod("con", construct!)],
        [MHMethod("li1", local_improve!)],
        Vector{MHMethod}();
        consider_initial_sol=false,
        log=true,
        checkit=true,
        titer=1000,
    )
    run!(alg)
    method_statistics(alg.scheduler)
    main_results(alg.scheduler)
    check(sol)

    write_solution("sol.txt", sol)
    return println(stderr, "spanning tree weight: $(sol.obj_val)")
    #println(stderr, "slack: $(resource_difference(inst, sol))")
end

#run_local_search()

function run_gvns()
    inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "100", "n100_k30_m2_s20019.txt"))
    #inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "200", "n200_k60_m3_s20021.txt"))
    #inst = RCkSTInstance(joinpath(@__DIR__, "..", "instances", "test", "1000", "n1000_k300_m3_s20021.txt"))
    #inst = RCkSTInstance(
    #    joinpath(
    #        @__DIR__, "..", "instances", "test", "10000", "n10000_k3000_m3_s20020.txt"
    #    ),
    #)

    sol::RCkSTSolution = RCkSTSolution(inst)
    seed!(123)

    # we apply a variable neighborhod search
    alg = GVNS(
        sol,
        [MHMethod("con", construct!)],
        [MHMethod("li1", local_improve!)],
        [MHMethod("sh1", shaking!, 1), MHMethod("sh2", shaking!, 2),
            MHMethod("sh3", shaking!, 3)];
        consider_initial_sol=false,
        log=true,
        #checkit=true,
        titer=1000,
    )
    run!(alg)
    method_statistics(alg.scheduler)
    main_results(alg.scheduler)
    check(sol)

    write_solution("sol.txt", sol)
    return println(stderr, "spanning tree weight: $(sol.obj_val)")
    #println(stderr, "slack: $(resource_difference(inst, sol))")
end

run_gvns()