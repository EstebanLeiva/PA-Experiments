using PulseAlgorithm
using DataStructures
using CSV
using Random
using DataFrames
using ProgressBars

function get_timeBudget(graph::Graph, start_node::Int, target_node::Int, α::Float64, γ::Float64, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
    shortest_mean_path = PA.dijkstra_between2Nodes(graph, start_node, target_node, "mean")
    shortest_cost_path = PA.dijkstra_between2Nodes(graph, start_node, target_node, "cost")
    
    cost_min_mean = 0.0
    for i in 1:length(shortest_mean_path)-1
        cost_min_mean = cost_min_mean + graph.nodes[shortest_mean_path[i]].links[shortest_mean_path[i+1]].cost 
    end

    cost_min_cost = 0.0
    for i in 1:length(shortest_cost_path)-1
        cost_min_cost = cost_min_cost + graph.nodes[shortest_cost_path[i]].links[shortest_cost_path[i+1]].cost 
    end
    
    mean, variance, covariance_term = PA.get_path_distribution(graph, shortest_mean_path, cov_dict)
    T_t_α = quantile(Normal(mean, √(variance+covariance_term)), α)

    mean, variance, covariance_term = PA.get_path_distribution(graph, shortest_cost_path, cov_dict)
    T_c_α = quantile(Normal(mean, √(variance+covariance_term)), α)

    T = T_t_α + (T_c_α - T_t_α) * (1 - γ) + 1e-3
    return T, shortest_mean_path, cost_min_mean, shortest_cost_path, cost_min_cost
end

function write_shortest_paths(graph::Graph, target_node::String, folder_path::String, network_name::String)
    target_node = graph.name_to_index[target_node]
    
    variance_costs = PA.dijkstra(graph, target_node, "variance")
    file_path = joinpath(folder_path, "variance_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(variance_costs = variance_costs), writeheader = false)

    mean_costs = PA.dijkstra(graph, target_node, "mean")
    file_path = joinpath(folder_path, "mean_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(mean_costs = mean_costs), writeheader = false)

    minimum_costs = PA.dijkstra(graph, target_node, "cost")
    file_path = joinpath(folder_path, "minimum_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(minimum_costs = minimum_costs), writeheader = false)
end

function preprocess_experiments(sp::SPulseGraph, folder_path::String, network_name::String, target_node::String)
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

    sp.source_node = rand(possible_start_nodes)[2].name     #set the source node randomly

    return sp.source_node
end

function run_experiments_time(graph::Graph, target_node::String, pulse::SPulseGraph, covariance_dict::DefaultDict, α::Float64, γ::Float64, folder_path::String, network_name::String, initial_bound::Bool)
    source_node = preprocess_experiments(pulse, folder_path, network_name, target_node) 
    T, shortest_mean_path, cost_min_mean, shortest_cost_path, cost_min_cost = get_timeBudget(graph, pulse.G.name_to_index[source_node], pulse.G.name_to_index[target_node], α, γ, covariance_dict)
    pulse.T_max = T
 
    mean_m, variance_m, covariance_term_m = PA.get_path_distribution(graph, shortest_mean_path, covariance_dict)
    reliability_shortest_mean_path_m = cdf(Normal(mean_m, √(variance_m+covariance_term_m)), T)

    mean_c, variance_c, covariance_term_c = PA.get_path_distribution(graph, shortest_cost_path, covariance_dict)
    reliability_shortest_cost_path_c = cdf(Normal(mean_c, √(variance_c+covariance_term_c)), T)

    if reliability_shortest_mean_path_m >= α && initial_bound
        elapsed_time = @elapsed begin
            PA.run_pulse(pulse, shortest_mean_path, cost_min_mean)
        end
    elseif reliability_shortest_cost_path_c >= α && initial_bound
        elapsed_time = @elapsed begin
            PA.run_pulse(pulse, shortest_cost_path, cost_min_cost)
        end
    else
        elapsed_time = @elapsed begin
            PA.run_pulse(pulse)
        end
    end
    return elapsed_time, pulse.instance_info, (source_node, target_node), pulse.T_max, pulse.optimal_path
end

function run_aggregated_experiments(graph::Graph, target_node::String, ρ::Float64, α::Float64, γ::Float64, max_depth::Int, folder_path::String, network_name::String, initial_bound::Bool, n::Int)
    covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)
    pulse = PA.create_SPulseGraph(graph, α, covariance_dict, target_node, target_node, 0.0)
    
    total_instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0, 
        "total_length_pruned_by_bounds" => 0,
        "total_length_pruned_by_feasibility" => 0,
        "number_nondominanted_paths" => 0, 
        "total_elapsed_time" => 0.0
    )
    
    for i in 1:n
        elapsed_time, instance_info, (start_node, target_node), T, optimal_path = run_experiments_time(graph, target_node, pulse, covariance_dict, α, γ, folder_path, network_name, initial_bound)
        total_instance_info["pruned_by_bounds"] += instance_info["pruned_by_bounds"]
        total_instance_info["pruned_by_feasibility"] += instance_info["pruned_by_feasibility"]
        total_instance_info["total_length_pruned_by_bounds"] += instance_info["total_length_pruned_by_bounds"]
        total_instance_info["total_length_pruned_by_feasibility"] += instance_info["total_length_pruned_by_feasibility"]
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += elapsed_time
    end
    return total_instance_info
end

function aggregate_experiment_results(sampled_keys::Vector{Int}, graph::Graph, ρ::Float64, α::Float64, γ::Float64, max_depth::Int, folder_path::String, network_name::String, n::Int)
    total_instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0, 
        "total_length_pruned_by_bounds" => 0,
        "total_length_pruned_by_feasibility" => 0,
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0
    )
    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info = run_aggregated_experiments(graph, graph.nodes[key].name, ρ, α, γ, max_depth, folder_path, network_name, true, n)
        total_instance_info["pruned_by_bounds"] += instance_info["pruned_by_bounds"]
        total_instance_info["pruned_by_feasibility"] += instance_info["pruned_by_feasibility"]
        total_instance_info["total_length_pruned_by_bounds"] += instance_info["total_length_pruned_by_bounds"]
        total_instance_info["total_length_pruned_by_feasibility"] += instance_info["total_length_pruned_by_feasibility"]
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += instance_info["total_elapsed_time"]
    end
    return total_instance_info
end