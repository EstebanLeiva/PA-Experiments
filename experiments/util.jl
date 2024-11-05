using ProgressBars

function get_timeBudget(graph::Graph, start_node::Int, target_node::Int, α::Float64, γ::Float64, cov_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64})
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

    mean, variance, covariance_term = PA.get_path_distribution(graph, shortest_mean_path, cov_dict)
    T_t_α = quantile(Normal(mean, √(variance + covariance_term)), α)

    mean, variance, covariance_term = PA.get_path_distribution(graph, shortest_cost_path, cov_dict)
    T_c_α = quantile(Normal(mean, √(variance + covariance_term)), α)

    T = T_t_α + (T_c_α - T_t_α) * (1 - γ) + 1e-3
    return T, shortest_mean_path, cost_min_mean, shortest_cost_path, cost_min_cost
end

function experiment_pasarp(graph::Graph, source_node::Int, target_node::Int, pulse::PaSarp, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, α::Float64, γ::Float64, initial_bound::Bool)
    T, shortest_mean_path, cost_min_mean, shortest_cost_path, cost_min_cost = get_timeBudget(graph, source_node, target_node, α, γ, covariance_dict)
    pulse.T_max = T

    mean_m, variance_m, covariance_term_m = PA.get_path_distribution(graph, shortest_mean_path, covariance_dict)
    reliability_shortest_mean_path_m = cdf(Normal(mean_m, √(variance_m + covariance_term_m)), T)

    mean_c, variance_c, covariance_term_c = PA.get_path_distribution(graph, shortest_cost_path, covariance_dict)
    reliability_shortest_cost_path_c = cdf(Normal(mean_c, √(variance_c + covariance_term_c)), T)

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
    return elapsed_time, pulse.instance_info
end

function n_experiments_pasarp(graph::Graph, target_node::Int, α::Float64, γ::Float64, cov_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, initial_bound::Bool, n::Int)
    pulse = PA.initialize_PaSarp(graph, α, cov_dict, graph.nodes[target_node].name, graph.nodes[target_node].name, 0.0)
    PA.preprocess!(pulse)
    total_instance_info = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "pruned_by_feasibility" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_feasibility" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}()
    )

    source_nodes = sample(collect(keys(graph.nodes)), n, replace=false)
    for source_node in source_nodes
        pulse.source_node = source_node
        elapsed_time, instance_info = experiment_pasarp(graph, source_node, target_node, pulse, cov_dict, α, γ, initial_bound)

        append!(total_instance_info["pruned_by_bounds"], instance_info["pruned_by_bounds"])
        append!(total_instance_info["pruned_by_feasibility"], instance_info["pruned_by_feasibility"])
        append!(total_instance_info["total_length_pruned_by_bounds"], instance_info["total_length_pruned_by_bounds"])
        append!(total_instance_info["total_length_pruned_by_feasibility"], instance_info["total_length_pruned_by_feasibility"])
        append!(total_instance_info["number_nondominanted_paths"], instance_info["number_nondominanted_paths"])
        append!(total_instance_info["total_elapsed_time"], elapsed_time)

    end
    return total_instance_info
end

function aggregate_experiments(sampled_keys::Vector{Int}, graph::Graph, α::Float64, γ::Float64, cov_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, n::Int)
    total_instance_info = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "pruned_by_feasibility" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_feasibility" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}()
    )
    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info = n_experiments_pasarp(graph, key, α, γ, cov_dict, true, n)

        total_instance_info["pruned_by_bounds"] = [total_instance_info["pruned_by_bounds"] ; instance_info["pruned_by_bounds"]]
        total_instance_info["pruned_by_feasibility"] = [total_instance_info["pruned_by_feasibility"] ; instance_info["pruned_by_feasibility"]]
        total_instance_info["total_length_pruned_by_bounds"] = [total_instance_info["total_length_pruned_by_bounds"] ; instance_info["total_length_pruned_by_bounds"]]
        total_instance_info["total_length_pruned_by_feasibility"] = [total_instance_info["total_length_pruned_by_feasibility"] ; instance_info["total_length_pruned_by_feasibility"]]
        total_instance_info["number_nondominanted_paths"] = [total_instance_info["number_nondominanted_paths"] ; instance_info["number_nondominanted_paths"]]
        total_instance_info["total_elapsed_time"] = [total_instance_info["total_elapsed_time"] ; instance_info["total_elapsed_time"]]

    end
    return total_instance_info
end