import MHLib
using Graphs: nv, ne, has_edge

include(joinpath(@__DIR__, "src", "RCkST.jl"))

# MHLib.check deliberately throws errors for invalid complete solutions.
# This helper ensures we catch the *expected* error, not an unrelated bug.
# Because it is a function, there are no Julia global soft-scope warnings.
function expect_error_contains(action::Function, expected_text::AbstractString)
    try
        action()
    catch err
        message = sprint(showerror, err)
        occursin(expected_text, message) || rethrow()
        println("   PASS: rejected with expected message: ", message)
        return nothing
    end
    error("Expected an error containing: $expected_text, but the call succeeded")
end

function check_instance_consistency(inst::RCkSTInstance)
    @assert size(inst.resources) == (inst.n, inst.m)
    @assert size(inst.edges) == (inst.edge_count, inst.m + 3)
    @assert length(inst.edges_hashed) == inst.edge_count
    @assert nv(inst.graph) == inst.n
    @assert ne(inst.graph) == inst.edge_count
    @assert all(r -> r >= 0, inst.resources)

    # Verify every edge ID, both orientations of dictionary lookup,
    # and Graphs.jl adjacency against the original edge matrix.
    for e in 1:inst.edge_count
        u, v = inst.edges[e, 1], inst.edges[e, 2]
        @assert 1 <= u <= inst.n && 1 <= v <= inst.n && u != v
        @assert edge_id(inst, u, v) == e
        @assert edge_id(inst, v, u) == e
        @assert has_edge(inst.graph, u, v)
        @assert all(c -> c >= 0, @view inst.edges[e, 3:end])
    end
    return nothing
end

function find_first_m2_file(folder::AbstractString)
    isdir(folder) || return nothing
    for (dir, _, files) in walkdir(folder)
        for file in sort(files)
            if endswith(file, ".txt") && occursin("_m2_", file)
                return joinpath(dir, file)
            end
        end
    end
    return nothing
end

function run_manual_checks()
    instance_dir = joinpath(@__DIR__, "instances")

    # ------------------------------------------------------------
    # 1. Original 6-node, 1-resource example: parser and graph
    # ------------------------------------------------------------
    println("\n1. Original example: parsing and lookup")
    original_path = joinpath(instance_dir, "example_instance.txt")
    @assert isfile(original_path) "Missing file: $original_path"
    inst = RCkSTInstance(original_path)
    @assert (inst.n, inst.edge_count, inst.k, inst.m) == (6, 8, 4, 1)
    @assert vec(inst.resources[:, 1]) == [3, 5, 2, 6, 1, 4]
    @assert vec(inst.edges[7, :]) == [4, 6, 1, 5]
    @assert edge_id(inst, 2, 4) == 4
    @assert edge_id(inst, 4, 2) == 4
    @assert isnothing(edge_id(inst, 1, 6))
    check_instance_consistency(inst)
    println("   PASS: parser, resources, graph, and dictionary")

    # ------------------------------------------------------------
    # 2. Valid complete solution: weight = 8, slack = [3]
    # ------------------------------------------------------------
    println("\n2. Original example: valid solution")
    s = RCkSTSolution(inst)
    set_edges!(s, [3, 4, 7])  # (2,3), (2,4), (4,6)
    @assert s.sel == inst.k - 1
    @assert selected_vertices(s) == Set([2, 3, 4, 6])
    @assert selected_is_connected(s, selected_vertices(s))
    @assert MHLib.obj(s) == 8
    @assert resource_slack(s) == [3]
    @assert MHLib.check(s) === nothing
    println("   PASS: objective=", MHLib.obj(s), ", slack=", resource_slack(s))

    # ------------------------------------------------------------
    # 3. Invalid complete solutions
    # ------------------------------------------------------------
    println("\n3. Invalid solutions must be rejected")
    bad_resources = set_edges!(RCkSTSolution(inst), [3, 5, 7])
    @assert MHLib.obj(bad_resources) == 6
    @assert resource_slack(bad_resources) == [-2]
    expect_error_contains("Resource constraints violated") do
        MHLib.check(bad_resources)
    end

    triangle = set_edges!(RCkSTSolution(inst), [3, 4, 5])
    @assert length(selected_vertices(triangle)) == 3
    expect_error_contains("Wrong number of selected vertices") do
        MHLib.check(triangle)
    end

    too_few = set_edges!(RCkSTSolution(inst), [3, 4])
    expect_error_contains("Wrong number of selected edges") do
        MHLib.check(too_few)
    end

    expect_error_contains("Duplicate edge IDs") do
        set_edges!(RCkSTSolution(inst), [3, 3, 7])
    end
    expect_error_contains("Invalid edge ID") do
        set_edges!(RCkSTSolution(inst), [3, 4, 999])
    end

    # This selection uses 5 distinct vertices and is disconnected.
    disconnected = set_edges!(RCkSTSolution(inst), [1, 7, 8])
    @assert !selected_is_connected(disconnected, selected_vertices(disconnected))
    println("   PASS: DFS detects disconnected selected edges")

    # With k=4 and 3 selected edges, it cannot simultaneously have
    # exactly k incident vertices and be disconnected in a simple graph.
    # Test MHLib.check's *connectivity branch* separately using k=5.
    mktempdir() do dir
        path = joinpath(dir, "connectivity_case.txt")
        write(path, """5 5 5 2
#node resources
10 10
10 10
10 10
10 10
10 10
#edges
1 2 1 1 1
2 3 1 1 1
1 3 1 1 1
4 5 1 1 1
3 4 1 1 1
""")
        five = RCkSTInstance(path)
        bad = set_edges!(RCkSTSolution(five), [1, 2, 3, 4])
        @assert length(selected_vertices(bad)) == five.k
        @assert bad.sel == five.k - 1
        @assert all(x -> x >= 0, resource_slack(bad))
        @assert !selected_is_connected(bad, selected_vertices(bad))
        expect_error_contains("Selected edges are disconnected") do
            MHLib.check(bad)
        end
    end

    # ------------------------------------------------------------
    # 4. Copy, copy!, empty!, cache invalidation
    # ------------------------------------------------------------
    println("\n4. Copy, reset, and cached objective")
    copy_s = copy(s)
    @assert copy_s !== s && copy_s.x !== s.x && copy_s.inst === s.inst
    empty!(copy_s)
    @assert copy_s.sel == 0 && s.sel == 3
    @assert MHLib.obj(copy_s) == 0
    copy!(copy_s, s)
    @assert MHLib.check(copy_s) === nothing
    set_edges!(copy_s, [3, 5, 7])
    @assert MHLib.obj(copy_s) == 6 && MHLib.obj(s) == 8
    println("   PASS: copies are independent and objective is updated")

    # ------------------------------------------------------------
    # 5. Solution file round-trip
    # ------------------------------------------------------------
    println("\n5. Save and reload a valid solution")
    mktempdir() do dir
        output = joinpath(dir, "solution.txt")
        write_solution(output, s)
        @assert readlines(output) == ["example_instance", "2 3", "2 4", "4 6"]
        loaded = read_solution(inst, output)
        @assert Set(loaded.x[1:loaded.sel]) == Set([3, 4, 7])
        @assert MHLib.obj(loaded) == 8
        @assert resource_slack(loaded) == [3]
        @assert MHLib.check(loaded) === nothing
        expect_error_contains("Resource constraints violated") do
            write_solution(joinpath(dir, "bad.txt"), bad_resources)
        end
    end
    println("   PASS: write_solution/read_solution agree")

    # ------------------------------------------------------------
    # 6. Small, explicitly constructed instance with TWO resources
    # ------------------------------------------------------------
    println("\n6. Small two-resource example: feasibility in both dimensions")
    small_path = joinpath(instance_dir, "example_instance_m2.txt")
    @assert isfile(small_path) "Save example_instance_m2.txt to: $small_path"
    two = RCkSTInstance(small_path)
    @assert (two.n, two.edge_count, two.k, two.m) == (5, 7, 3, 2)
    @assert two.resources == [2 7; 4 2; 1 5; 5 1; 2 2]
    check_instance_consistency(two)

    feasible_two = set_edges!(RCkSTSolution(two), [1, 2])
    @assert selected_vertices(feasible_two) == Set([1, 2, 3])
    @assert MHLib.obj(feasible_two) == 6
    @assert resource_slack(feasible_two) == [2, 8]
    @assert MHLib.check(feasible_two) === nothing
    println("   PASS: feasible m=2 tree: weight=6, slack=[2, 8]")

    # Both selections are connected trees with k=3 vertices,
    # but one violates only resource 1 and one violates only resource 2.
    bad_first = set_edges!(RCkSTSolution(two), [1, 3])
    @assert resource_slack(bad_first) == [-3, 3]
    expect_error_contains("Resource constraints violated") do
        MHLib.check(bad_first)
    end
    bad_second = set_edges!(RCkSTSolution(two), [2, 4])
    @assert resource_slack(bad_second) == [7, -1]
    expect_error_contains("Resource constraints violated") do
        MHLib.check(bad_second)
    end
    println("   PASS: both separate resource violations are caught")

    # ------------------------------------------------------------
    # 7. REAL course instance with m=2 (screenshot: test/100)
    # ------------------------------------------------------------
    println("\n7. Real test instance (100 vertices, 2 resources)")
    real_path = joinpath(instance_dir, "test", "100", "n100_k30_m2_s20001.txt")
    @assert isfile(real_path) "Expected course instance was not found: $real_path"
    real = RCkSTInstance(real_path)
    @assert (real.n, real.k, real.m) == (100, 30, 2)
    check_instance_consistency(real)
    println("   PASS: n=", real.n, ", edges=", real.edge_count,
            ", k=", real.k, ", m=", real.m)
    println("   NOTE: This verifies parsing and graph consistency, not a feasible tree;")
    println("         no known solution was provided for the real instance.")

    # ------------------------------------------------------------
    # 8. REAL training instance, if an m=2 file exists
    # ------------------------------------------------------------
    println("\n8. Real training instance with two resources")
    train_path = find_first_m2_file(joinpath(instance_dir, "train"))
    if isnothing(train_path)
        println("   SKIPPED: no *_m2_*.txt file found below instances/train/")
    else
        train = RCkSTInstance(train_path)
        @assert train.m == 2
        check_instance_consistency(train)
        println("   PASS: ", basename(train_path), " (n=", train.n,
                ", edges=", train.edge_count, ", k=", train.k, ")")
    end

    println("\nALL EXECUTED MANUAL CHECKS PASSED")
    return nothing
end

run_manual_checks()
