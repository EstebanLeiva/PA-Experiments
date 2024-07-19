using Distributions
include("util.jl")

function modified_dfs(graph::Graph, start_link::Tuple{Int, Int}, max_depth::Int, depth::Int, visited_pairlinks::Dict{Tuple{Int, Int}, Int}, previous_node::Int)
    if depth > max_depth
        return visited_pairlinks
    end

    for adjacent in keys(graph.nodes[start_link[2]].links)
        if previous_node != adjacent
            if haskey(visited_pairlinks, (start_link[2], adjacent))
                visited_pairlinks[(start_link[2], adjacent)] = min(visited_pairlinks[(start_link[2], adjacent)], depth)
            else
                visited_pairlinks[(start_link[2], adjacent)] = depth
            end
            modified_dfs(graph, (start_link[2], adjacent), max_depth, depth + 1, visited_pairlinks, adjacent)
        end
    end       

    return visited_pairlinks
end

function get_neg_covariance_dict(graph::Graph, ρ::Float64, max_depth::Int)
    covariance_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0) #default value of dic is 0.0
    links = PA.get_links_info(graph)
    for link in keys(links)
        visited_pairlinks = modified_dfs(graph, link, max_depth, 1, Dict{Tuple{Int, Int}, Int}(), -1)
        for pairlink in keys(visited_pairlinks)
            if pairlink != link 
                #do a coin toss with 30% probability of being negative
                if rand() < 0.7
                    covariance_dict[(link[1], link[2], pairlink[1], pairlink[2])] = (ρ/(visited_pairlinks[pairlink])) * √(links[(pairlink)][3]) * √(links[link][3]) #Corredor structure of ρ's
                else
                    covariance_dict[(link[1], link[2], pairlink[1], pairlink[2])] = - (ρ/(visited_pairlinks[pairlink])) * √(links[(pairlink)][3]) * √(links[link][3]) #Corredor structure of ρ's
                end
            end
            if pairlink == link
                covariance_dict[(link[1],link[2],pairlink[1],pairlink[2])] = 1.0
            end
        end
    end
    return covariance_dict
end

function run_neg_aggregated_experiments(graph::Graph, target_node::String, ρ::Float64, α::Float64, γ::Float64, max_depth::Int, folder_path::String, network_name::String, initial_bound::Bool, n::Int)
    covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)
    pulse = PA.initialize(graph, α, covariance_dict, target_node, target_node, 0.0)
    
    total_instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0, 
        "total_length_pruned_by_bounds" => 0,
        "total_length_pruned_by_feasibility" => 0,
        "number_nondominanted_paths" => 0, 
        "total_elapsed_time" => 0.0, 
        "post_reliability" => 0.0
    )
    
    for i in 1:n
        elapsed_time, instance_info, (start_node, _), T, optimal_path = run_experiments_time(graph, target_node, pulse, covariance_dict, α, γ, folder_path, network_name, initial_bound)
        mean, variance, cov = PA.get_path_distribution(graph, optimal_path, covariance_dict)
        dist = Normal(mean, √(variance + cov))
        reliability = cdf(dist, T)
        total_instance_info["pruned_by_bounds"] += instance_info["pruned_by_bounds"]
        total_instance_info["pruned_by_feasibility"] += instance_info["pruned_by_feasibility"]
        total_instance_info["total_length_pruned_by_bounds"] += instance_info["total_length_pruned_by_bounds"]
        total_instance_info["total_length_pruned_by_feasibility"] += instance_info["total_length_pruned_by_feasibility"]
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += elapsed_time
        total_instance_info["post_reliability"] += reliability
    end
    return total_instance_info
end

function aggregate_neg_experiment_results(sampled_keys::Vector{Int}, graph::Graph, ρ::Float64, α::Float64, γ::Float64, max_depth::Int, folder_path::String, network_name::String, n::Int)
    total_instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0, 
        "total_length_pruned_by_bounds" => 0,
        "total_length_pruned_by_feasibility" => 0,
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0, 
        "post_reliability" => 0.0
    )
    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info = run_neg_aggregated_experiments(graph, graph.nodes[key].name, ρ, α, γ, max_depth, folder_path, network_name, true, n)
        total_instance_info["pruned_by_bounds"] += instance_info["pruned_by_bounds"]
        total_instance_info["pruned_by_feasibility"] += instance_info["pruned_by_feasibility"]
        total_instance_info["total_length_pruned_by_bounds"] += instance_info["total_length_pruned_by_bounds"]
        total_instance_info["total_length_pruned_by_feasibility"] += instance_info["total_length_pruned_by_feasibility"]
        total_instance_info["number_nondominanted_paths"] += instance_info["number_nondominanted_paths"]
        total_instance_info["total_elapsed_time"] += instance_info["total_elapsed_time"]
        total_instance_info["post_reliability"] += instance_info["post_reliability"]
    end
    return total_instance_info
end