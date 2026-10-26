
import MHLib

# Load implementation
include(joinpath(@__DIR__, "src", "RCkST.jl"))

# Load the example instance
inst = RCkSTInstance(
    joinpath(@__DIR__, "instances", "example_instance.txt")
)

# ========================================
# TEST 1: INSTANCE PARSING
# ========================================

println("TEST 1: Reading the instance")

println("Vertices: ", inst.n)
println("Edges: ", inst.edge_count)
println("k: ", inst.k)
println("Resources: ", inst.m)

@assert inst.n == 6
@assert inst.edge_count == 8
@assert inst.k == 4
@assert inst.m == 1

println("PASS: Instance parsing")


# ========================================
# TEST 2: VALID SOLUTION
# ========================================

println("\nTEST 2: Valid solution")

s = RCkSTSolution(inst)

# Select edges:
# 3: (2,3), weight 2, consumption 5
# 4: (2,4), weight 5, consumption 4
# 7: (4,6), weight 1, consumption 5

set_edges!(s, [3, 4, 7])

println("Selected edges: ", s.x[1:s.sel])
println("Selected vertices: ", selected_vertices(s))
println("Objective: ", MHLib.obj(s))
println("Resource slack: ", resource_slack(s))

@assert selected_vertices(s) == Set([2, 3, 4, 6])
@assert MHLib.obj(s) == 8
@assert resource_slack(s) == [3]

# Check complete solution
MHLib.check(s)

println("PASS: Solution is feasible")



# ========================================
# TEST 3: INFEASIBLE SOLUTION
# ========================================

println("\nTEST 3: Infeasible solution")

let
    bad = RCkSTSolution(inst)

    # Cheaper tree, but consumes too many resources
    set_edges!(bad, [3, 5, 7])

    println("Objective: ", MHLib.obj(bad))
    println("Resource slack: ", resource_slack(bad))

    @assert MHLib.obj(bad) == 6
    @assert resource_slack(bad) == [-2]

    rejected = false

    try
        MHLib.check(bad)

    catch err
        message = sprint(showerror, err)

        if occursin("Resource constraints violated", message)
            rejected = true
            println("Expected error: ", message)
        else
            rethrow()
        end
    end

    @assert rejected "Invalid solution was accepted"

    println("PASS: Invalid solution rejected")
end

println("\nALL TESTS PASSED")
