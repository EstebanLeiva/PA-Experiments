using PulseAlgorithm
using DataStructures
using CSV
using Random
using DataFrames
using ProgressBars

function get_initial_paths(graph::Graph, start_node::Int, target_node::Int, α::Float64, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
    shortest_mean_path = PA.dijkstra_between_nodes(graph, start_node, target_node, "mean")
    shortest_cost_path = PA.dijkstra_between_nodes(graph, start_node, target_node, "cost")
    
    cost_min_mean = 0.0
    for i in 1:length(shortest_mean_path)-1
        cost_min_mean = cost_min_mean + graph.nodes[shortest_mean_path[i]].links[shortest_mean_path[i+1]].cost 
    end

    cost_min_cost = 0.0
    for i in 1:length(shortest_cost_path)-1
        cost_min_cost = cost_min_cost + graph.nodes[shortest_cost_path[i]].links[shortest_cost_path[i+1]].cost 
    end

    return shortest_mean_path, cost_min_mean, shortest_cost_path, cost_min_cost
end

function preprocess_experiments(sp::SdrsppPulse, folder_path::String, network_name::String, target_node::String)
    file_path = joinpath(folder_path, "minimum_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.minimum_costs =  data[:, 1] |> collect

    file_path = joinpath(folder_path, "mean_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.mean_costs =  data[:, 1] |> collect

    file_path = joinpath(folder_path, "variance_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.variance_costs =  data[:, 1] |> collect

    nodes = sort(collect(sp.G.nodes), by=x->sp.minimum_costs[x[1]], rev = true)
    possible_start_nodes = filter(x->sp.minimum_costs[x[1]] != Inf, nodes)

    sp.source_node = sp.G.name_to_index[rand(possible_start_nodes)[2].name]  #set the source node randomly

    return sp.source_node
end

function run_experiments_time(graph::Graph, target_node::String, pulse::SdrsppPulse, covariance_dict::DefaultDict, α::Float64, folder_path::String, network_name::String, initial_bound::Bool)
    source_node = preprocess_experiments(pulse, folder_path, network_name, target_node) 
    target_node = graph.name_to_index[target_node]
    shortest_mean_path, _, _, _ = get_initial_paths(graph, source_node, target_node, α, covariance_dict)
 
    mean_m, variance_m, covariance_term_m = PA.get_path_distribution(graph, shortest_mean_path, covariance_dict)
    quantile_m = quantile(Normal(mean_m, √(variance_m+covariance_term_m)), α)

    if initial_bound
        elapsed_time = @elapsed begin
            PA.run_pulse(pulse, shortest_mean_path, quantile_m)
        end
    end

    return elapsed_time, pulse.instance_info, (source_node, target_node), pulse.optimal_path
end

function run_aggregated_experiments(graph::Graph, target_node::String, ρ::Float64, α::Float64, max_depth::Int, folder_path::String, network_name::String, initial_bound::Bool, n::Int)
    covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)
    pulse = PA.initialize(graph, α, covariance_dict, target_node, target_node)
    
    total_instance_info = Dict(
        "pruned_by_bounds" => 0,
        "total_length_pruned_by_bounds" => 0,
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0
        )
    
    for i in 1:n
        elapsed_time, instance_info, (start_node, _), optimal_path = run_experiments_time(graph, target_node, pulse, covariance_dict, α, folder_path, network_name, initial_bound)
        total_instance_info["pruned_by_bounds"] += instance_info["pruned_by_bounds"]
        total_instance_info["total_length_pruned_by_bounds"] += instance_info["total_length_pruned_by_bounds"]
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += elapsed_time
    end
    return total_instance_info
end

function aggregate_experiment_results(sampled_keys::Vector{Int}, graph::Graph, ρ::Float64, α::Float64, max_depth::Int, folder_path::String, network_name::String, n::Int)
    total_instance_info = Dict(
        "pruned_by_bounds" => 0,
        "total_length_pruned_by_bounds" => 0,
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0
    )
    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info = run_aggregated_experiments(graph, graph.nodes[key].name, ρ, α, max_depth, folder_path, network_name, true, n)
        total_instance_info["pruned_by_bounds"] += instance_info["pruned_by_bounds"]
        total_instance_info["total_length_pruned_by_bounds"] += instance_info["total_length_pruned_by_bounds"]
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += instance_info["total_elapsed_time"]
    end
    return total_instance_info
end

### Erspa ###
function erspa_preprocess_experiments(graph::Graph, covariance_dict::DefaultDict, folder_path::String, network_name::String, target_node::String, α::Float64)
    file_path = joinpath(folder_path, "minimum_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    minimum_costs =  data[:, 1] |> collect

    nodes = sort(collect(graph.nodes), by=x->minimum_costs[x[1]], rev = true)
    possible_start_nodes = filter(x->minimum_costs[x[1]] != Inf, nodes)

    source_node = rand(possible_start_nodes)[2].name 

    folder_path = folder_path[1:end-14]
    file_path = joinpath(folder_path, "node_coordinates_" * network_name * ".tntp")
    node_coordinates = Vector{Tuple{Float64, Float64}}()
    open(file_path, "r") do file
        readline(file)
        for line in eachline(file)
            parts = split(line, '\t')
            x = parse(Float64, parts[2])
            y = parse(Float64, parts[3])
            push!(node_coordinates, (x, y))
        end
    end

    erspa = PA.initialize(graph, α, covariance_dict, node_coordinates, source_node, target_node)
    PA.preprocess!(erspa)
    
    return erspa
end

function erspa_run_experiments_time(erspa::ErspaStar)
    elapsed_time = @elapsed begin
        PA.run_erspa(erspa)
    end
    return elapsed_time, erspa.instance_info, (erspa.source_node, erspa.target_node), erspa.optimal_path
end

function erspa_run_aggregated_experiments(graph::Graph, target_node::String, ρ::Float64, α::Float64, max_depth::Int, folder_path::String, network_name::String, n::Int)
    covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)
    
    total_instance_info = Dict(
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0
        )
    
    for i in 1:n
        erspa = erspa_preprocess_experiments(graph, covariance_dict, folder_path, network_name, target_node, α) 
        elapsed_time, instance_info, (start_node, _), optimal_path = erspa_run_experiments_time(erspa)
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += elapsed_time
    end
    return total_instance_info
end

function erspa_aggregate_experiment_results(sampled_keys::Vector{Int}, graph::Graph, ρ::Float64, α::Float64, max_depth::Int, folder_path::String, network_name::String, n::Int)
    total_instance_info = Dict(
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0
    )
    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info = erspa_run_aggregated_experiments(graph, graph.nodes[key].name, ρ, α, max_depth, folder_path, network_name, n)
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += instance_info["total_elapsed_time"]
    end
    return total_instance_info
end