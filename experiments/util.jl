using ProgressBars

function get_path_distribution(graph::PA.Graph, path::Vector{Int}, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
    mean = 0.0
    variance = 0.0
    covariance_term = 0.0
    for i in 1:length(path)-1
        mean += graph.nodes[path[i]].links[path[i+1]].random["time"]["mean"]
        variance += graph.nodes[path[i]].links[path[i+1]].random["time"]["variance"]
        for ii in i + 1:length(path)-1
            covariance_term += 2*cov_dict[(path[i], path[i+1], path[ii], path[ii+1])]
        end
    end
    return mean, variance, covariance_term
end

function get_covariance_term(covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, reachable_node::Int, path::Vector{Int})
    n = length(path)
    if n > 1
        last_node = path[end]
        current_node = reachable_node
        covariance_sum = 0.0
        for i in 1:n-1
            covariance_sum += 2 * covariance_dict[(path[i], path[i+1], last_node, current_node)]
        end
        return covariance_sum
    else
        return 0.0
    end
end

function get_timeBudget(graph::PA.Graph, start_node::Int, target_node::Int, α::Float64, γ::Float64, cov_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64})
    shortest_mean_path = PA.dijkstra(graph, start_node, target_node, "time", "mean")
    shortest_cost_path = PA.dijkstra(graph, start_node, target_node, "cost")

    cost_min_mean = 0.0
    for i in 1:length(shortest_mean_path)-1
        cost_min_mean = cost_min_mean + graph.nodes[shortest_mean_path[i]].links[shortest_mean_path[i+1]].deterministic["cost"]
    end

    cost_min_cost = 0.0
    for i in 1:length(shortest_cost_path)-1
        cost_min_cost = cost_min_cost + graph.nodes[shortest_cost_path[i]].links[shortest_cost_path[i+1]].deterministic["cost"]
    end

    mean, variance, covariance_term = get_path_distribution(graph, shortest_mean_path, cov_dict)
    T_t_α = quantile(Normal(mean, √(variance + covariance_term)), α)

    mean, variance, covariance_term = get_path_distribution(graph, shortest_cost_path, cov_dict)
    T_c_α = quantile(Normal(mean, √(variance + covariance_term)), α)

    T = T_t_α + (T_c_α - T_t_α) * (1 - γ) + 1e-3
    return T, shortest_mean_path, cost_min_mean, shortest_cost_path, cost_min_cost
end

function experiment_PaSarp(pulse::PA.Pulse, info_update::Function, pruning_functions::Vector{Function}, pulse_score::Function, initial_bound::Bool)
    preprocessing_time = @elapsed begin
        preprocess!(pulse)
    end
    graph = pulse.problem.graph
    source_node = pulse.problem.source_node
    target_node = pulse.problem.target_node
    α = pulse.problem.constants["alpha"]
    γ = pulse.problem.constants["gamma"]
    covariance_dict = pulse.problem.graph.covariance["time"]
    T, shortest_mean_path, cost_min_mean, shortest_cost_path, cost_min_cost = get_timeBudget(graph, source_node, target_node, α, γ, covariance_dict)
    pulse.problem.constants["T_max"] = T

    mean_m, variance_m, covariance_term_m = get_path_distribution(graph, shortest_mean_path, covariance_dict)
    reliability_shortest_mean_path_m = cdf(Normal(mean_m, √(variance_m + covariance_term_m)), T)

    mean_c, variance_c, covariance_term_c = get_path_distribution(graph, shortest_cost_path, covariance_dict)
    reliability_shortest_cost_path_c = cdf(Normal(mean_c, √(variance_c + covariance_term_c)), T)

    if reliability_shortest_mean_path_m >= α && initial_bound
        elapsed_time = @elapsed begin
            PA.run_pulse!(pulse, info_update, pruning_functions, pulse_score, shortest_mean_path, cost_min_mean)
        end
    elseif reliability_shortest_cost_path_c >= α && initial_bound
        elapsed_time = @elapsed begin
            PA.run_pulse!(pulse, info_update, pruning_functions, pulse_score, shortest_cost_path, cost_min_cost)
        end
    else
        elapsed_time = @elapsed begin
            PA.run_pulse!(pulse, info_update, pruning_functions, pulse_score)
        end
    end
    return elapsed_time, preprocessing_time, pulse.instance_info
end

function n_experiments_PaSarp(graph::PA.Graph, target_node::Int, parameters::PA.Parameters, constants::Dict{String, Float64}, pruning_functions::Vector{Function}, info_update::Function, pulse_score::Function, initial_bound::Bool, n::Int)
    total_instance_info = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "pruned_by_feasibility" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_feasibility" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(),
        "total_preprocessing_time" => Vector{Float64}()
    )

    source_nodes = sample(collect(keys(graph.nodes)), n, replace=false)
    for source_node in source_nodes
        problem = PA.Problem(graph, source_node, target_node, true, constants)
        pulse = PA.Pulse(problem, parameters)
        # Initialize instance_info
        pulse.instance_info["pruned_by_bounds"] = 0
        pulse.instance_info["pruned_by_feasibility"] = 0
        pulse.instance_info["total_length_pruned_by_bounds"] = 0
        pulse.instance_info["total_length_pruned_by_feasibility"] = 0
        pulse.instance_info["number_nondominanted_paths"] = 0

        # These are the error values, that will not be taken into account for the analysis
        elapsed_time = 123456789
        preprocessing_time = 123456789
        instance_info = Dict{String, Float64}()
        instance_info["pruned_by_bounds"] = 123456789
        instance_info["pruned_by_feasibility"] = 123456789
        instance_info["total_length_pruned_by_bounds"] = 123456789
        instance_info["total_length_pruned_by_feasibility"] = 123456789
        instance_info["number_nondominanted_paths"] = 123456789

        try
            elapsed_time, preprocessing_time, instance_info = experiment_PaSarp(pulse, info_update, pruning_functions, pulse_score, initial_bound)
        catch e
            @warn "Non-reachable node" exception=(e, catch_backtrace())
        end

        append!(total_instance_info["pruned_by_bounds"], instance_info["pruned_by_bounds"])
        append!(total_instance_info["pruned_by_feasibility"], instance_info["pruned_by_feasibility"])
        append!(total_instance_info["total_length_pruned_by_bounds"], instance_info["total_length_pruned_by_bounds"])
        append!(total_instance_info["total_length_pruned_by_feasibility"], instance_info["total_length_pruned_by_feasibility"])
        append!(total_instance_info["number_nondominanted_paths"], instance_info["number_nondominanted_paths"])
        append!(total_instance_info["total_elapsed_time"], elapsed_time)
        append!(total_instance_info["total_preprocessing_time"], preprocessing_time)

    end
    return total_instance_info
end

function aggregate_experiments(sampled_keys::Vector{Int}, graph::PA.Graph, parameters::PA.Parameters, constants::Dict{String, Float64}, pruning_functions::Vector{Function}, info_update::Function, pulse_score::Function, n::Int)
    total_instance_info = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "pruned_by_feasibility" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_feasibility" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(),
        "total_preprocessing_time" => Vector{Float64}()
    )
    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info = n_experiments_PaSarp(graph, key, parameters, constants, pruning_functions, info_update, pulse_score, true, n)

        total_instance_info["pruned_by_bounds"] = [total_instance_info["pruned_by_bounds"] ; instance_info["pruned_by_bounds"]]
        total_instance_info["pruned_by_feasibility"] = [total_instance_info["pruned_by_feasibility"] ; instance_info["pruned_by_feasibility"]]
        total_instance_info["total_length_pruned_by_bounds"] = [total_instance_info["total_length_pruned_by_bounds"] ; instance_info["total_length_pruned_by_bounds"]]
        total_instance_info["total_length_pruned_by_feasibility"] = [total_instance_info["total_length_pruned_by_feasibility"] ; instance_info["total_length_pruned_by_feasibility"]]
        total_instance_info["number_nondominanted_paths"] = [total_instance_info["number_nondominanted_paths"] ; instance_info["number_nondominanted_paths"]]
        total_instance_info["total_elapsed_time"] = [total_instance_info["total_elapsed_time"] ; instance_info["total_elapsed_time"]]
        total_instance_info["total_preprocessing_time"] = [total_instance_info["total_preprocessing_time"] ; instance_info["total_preprocessing_time"]]

    end
    return total_instance_info
end

function modified_dfs(graph::PA.Graph, start_link::Tuple{Int, Int}, max_depth::Int, depth::Int, visited_pairlinks::Dict{Tuple{Int, Int}, Int}, previous_node::Int)
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

function get_covariance_dict(graph::PA.Graph, random_variable::String, ρ::Float64, max_depth::Int)
    covariance_dict = Dict{String, DefaultDict{Tuple{Int, Int, Int, Int}, Float64}}() 
    covariance_dict[random_variable] = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0)
    links = PA.get_links_info(graph)
    for link in keys(links)
        visited_pairlinks = modified_dfs(graph, link, max_depth, 1, Dict{Tuple{Int, Int}, Int}(), -1)
        for pairlink in keys(visited_pairlinks)
            if pairlink != link 
                covariance_dict[random_variable][(link[1], link[2], pairlink[1], pairlink[2])] = (ρ/(visited_pairlinks[pairlink])) * √(links[(pairlink)][2][random_variable]["variance"]) * √(links[link][2][random_variable]["variance"]) #Corredor structure of ρ's
            end
            if pairlink == link
                covariance_dict[random_variable][(link[1],link[2],pairlink[1],pairlink[2])] = 1.0
            end
        end
    end
    return covariance_dict
end