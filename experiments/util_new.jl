using ProgressBars
using PulseAlgorithm
const PA = PulseAlgorithm

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

function sample_od_pairs(graph::PA.Graph, 
                            num_target_nodes::Int,
                            num_source_nodes::Int)

    sampled_keys = sample(collect(keys(graph.nodes)), num_target_nodes, replace=false)
    pair_dict = Dict{Int, Vector{Int}}()
    for key in sampled_keys
        pair_dict[key] = sample(collect(keys(graph.nodes)), num_source_nodes, replace=false)
    end
    return pair_dict
end

function experiments_sarp(pair_dict::Dict{Int, Vector{Int}},
                          graph::PA.Graph, 
                          parameters:: PA.Parameters,
                          params_weight::PA.Parameters,
                          constants::Dict{String, Float64},
                          pruning_functions,
                          info_update,
                          pulse_score,
                          file_name)

    total_instance_info = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "pruned_by_feasibility" => Vector{Int}(),
        "length_pruned_by_bounds" => Vector{Int}(),
        "length_pruned_by_feasibility" => Vector{Int}(),
        "number_nondominated_paths" => Vector{Int}(),
        "elapsed_time" => Vector{Float64}(),
        "preprocessing_time" => Vector{Float64}(),
        "optimal_path_length" => Vector{Int}(),
        "distance_source_target" => Vector{Int}(),
        "distance" => Vector{Int}(),
        "timeout" => Vector{Bool}(),
        "optimal_path" => Vector{String}(),
        "target" => Vector{Int}(),
        "source" => Vector{Int}()
    )

    for target in ProgressBar(keys(pair_dict))
        for source in pair_dict[target]
            problem = PA.Problem(graph, source, target, true, constants)
            pulse_weight = PA.Pulse(problem, params_weight)
            pulse = PA.Pulse(problem, parameters)
            pulse.parameters.timeout = false
            # Initialize instance_info
            pulse.instance_info["pruned_by_bounds"] = 0
            pulse.instance_info["pruned_by_feasibility"] = 0
            pulse.instance_info["length_pruned_by_bounds"] = 0
            pulse.instance_info["length_pruned_by_feasibility"] = 0
            pulse.instance_info["number_nondominated_paths"] = 0
            
            try

                elapsed_time, preprocessing_time, optimal_path, optimal_path_length, distance, instance_info, timeout = one_sarp(pulse, 
                                                                                                pulse_weight,
                                                                                                info_update, 
                                                                                                pruning_functions, 
                                                                                                pulse_score
                                                                                                )
                append!(total_instance_info["pruned_by_bounds"], instance_info["pruned_by_bounds"])
                append!(total_instance_info["pruned_by_feasibility"], instance_info["pruned_by_feasibility"])
                append!(total_instance_info["length_pruned_by_bounds"], instance_info["length_pruned_by_bounds"])
                append!(total_instance_info["length_pruned_by_feasibility"], instance_info["length_pruned_by_feasibility"])
                append!(total_instance_info["number_nondominated_paths"], instance_info["number_nondominated_paths"])
                append!(total_instance_info["elapsed_time"], elapsed_time)
                append!(total_instance_info["preprocessing_time"], preprocessing_time)
                append!(total_instance_info["optimal_path_length"], optimal_path_length)
                append!(total_instance_info["distance"], distance)
                append!(total_instance_info["timeout"], timeout)
                append!(total_instance_info["target"], target)
                append!(total_instance_info["source"], source)

            catch e
                @warn "Non-reachable node" exception=(e, catch_backtrace())
            end
        end
    end

    length_pruned_by_bounds = (total_instance_info["length_pruned_by_bounds"] ./ total_instance_info["pruned_by_bounds"])
    length_pruned_by_feasibility = (total_instance_info["length_pruned_by_feasibility"] ./ total_instance_info["pruned_by_feasibility"])

    df_prune = DataFrame(
        length_pruned_by_bounds = length_pruned_by_bounds,
        length_pruned_by_feasibility = length_pruned_by_feasibility
    )

    #println(total_instance_info["optimal_path"])
    df = DataFrame(
        elapsed_time = total_instance_info["elapsed_time"],
        preprocessing_time = total_instance_info["preprocessing_time"],
        optimal_path_length = total_instance_info["optimal_path_length"],
        number_nondominated_paths = total_instance_info["number_nondominated_paths"],
        pruned_by_bounds = total_instance_info["pruned_by_bounds"],
        pruned_by_feasibility = total_instance_info["pruned_by_feasibility"],
        distance = total_instance_info["distance"],
        timeout = total_instance_info["timeout"],
        target = total_instance_info["target"],
        source = total_instance_info["source"]
    )
    name = file_name * "_pruning.csv"
    CSV.write(name, df_prune)
    name = file_name * "_experiment.csv"
    CSV.write(name, df)
end

function one_sarp(pulse, pulse_weight, info_update, pruning_functions, pulse_score)
    preprocess!(pulse_weight)
    distance = pulse_weight.preprocessing.deterministic["weight"][pulse.problem.source_node]

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
    pulse.parameters.max_pulse_depth = length(shortest_mean_path)

    mean_m, variance_m, covariance_term_m = get_path_distribution(graph, shortest_mean_path, covariance_dict)
    reliability_shortest_mean_path_m = cdf(Normal(mean_m, √(variance_m + covariance_term_m)), T)

    mean_c, variance_c, covariance_term_c = get_path_distribution(graph, shortest_cost_path, covariance_dict)
    reliability_shortest_cost_path_c = cdf(Normal(mean_c, √(variance_c + covariance_term_c)), T)

    
    if reliability_shortest_mean_path_m >= α
        elapsed_time = @elapsed begin
            run_pulse!(pulse, info_update, pruning_functions, pulse_score, init_optimal_path=shortest_mean_path, init_objective=cost_min_mean)
        end
    elseif reliability_shortest_cost_path_c >= α
        elapsed_time = @elapsed begin
            run_pulse!(pulse, info_update, pruning_functions, pulse_score, init_optimal_path=shortest_cost_path, init_objective=cost_min_cost)
        end
    else
        elapsed_time = @elapsed begin
            run_pulse!(pulse, info_update, pruning_functions, pulse_score)
        end
    end

    return elapsed_time, preprocessing_time, pulse.optimal_path, length(pulse.optimal_path), distance, pulse.instance_info, pulse.parameters.timeout
end