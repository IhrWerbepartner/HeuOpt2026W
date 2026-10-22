using MHLib

export RCkSTInstance, RCkSTSolution

"""
- `n`: number of vertices
- `m`: number of resources per vertex
- `k`: edges in spanning tree
- `resoucres`: resources available at each vertex resources[v]
- `edges`: edges in the graph with corresponding costs [1..m]
"""
struct RCkSTInstance
    file_name::String
    n::Int
    m::Int
    k::Int
    resources::Matrix{Int}
    edges::Matrix{Int}
end

struct Edge
    v0::Int
    v1::Int
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
        error("invalid number of values in RCkST instance: $(length(all_values)), expected: $expected_values")
    end

    resources_end_index::Int = 5 + n * m - 1
    resources = reshape(all_values[5:resources_end_index], (n, m))
    edges = reshape(all_values[(resources_end_index+1):(resources_end_index+edge_count*(3+m))],
        (edge_count, 3 + m))

    # TODO: figure out how to rescale the resources but avoid floating point numbers (possibly (?)) for performance reasons
    RCkSTInstance(file_name, n, m, k, resources, edges)
end

# test reading instances is working
RCkSTInstance("../instances/example_instance.txt")


"""
TODO
"""
function validate_solution(solution::RCkSTInstance)
end