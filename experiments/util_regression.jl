using ProgressBars

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
    pulse.parameters.max_pulse_depth = length(shortest_mean_path)

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
    distance = length(PA.dijkstra(graph, source_node, target_node, "weight"))
    return elapsed_time, preprocessing_time, T, distance
end

function n_experiments_PaSarp(graph::PA.Graph, 
                            target_node::Int, 
                            parameters::PA.Parameters, 
                            constants::Dict{String, Real}, 
                            pruning_functions::Vector{Function}, 
                            info_update::Function, 
                            pulse_score::Function, 
                            initial_bound::Bool, 
                            n::Int)

    info = Dict(
        "elapsed_time" => Float64[],
        "preprocessing_time" => Float64[],
        "time_budget" => Float64[], 
        "alpha" => Float64[],
        "gamma" => Float64[],
        "rho" => Float64[],
        "max_depth" => Int[],
        "distance" => Int[],
    )
    source_nodes = sample(collect(keys(graph.nodes)), n, replace=false)
    for source_node in source_nodes
        problem = PA.Problem(graph, source_node, target_node, true, constants)
        pulse = PA.Pulse(problem, parameters)
    
        # These are the error values, that will not be taken into account for the analysis
        elapsed_time = 123456789
        preprocessing_time = 123456789

        try
            elapsed_time, preprocessing_time, time_budget, distance = experiment_PaSarp(pulse, info_update, pruning_functions, pulse_score, initial_bound)
            append!(info["elapsed_time"], elapsed_time)
            append!(info["preprocessing_time"], preprocessing_time)
            append!(info["time_budget"], time_budget)
            append!(info["distance"], distance)
            append!(info["alpha"], constants["alpha"])
            append!(info["gamma"], constants["gamma"])
            append!(info["rho"], constants["rho"])
            append!(info["max_depth"], constants["max_depth"])
        catch e
            @warn "Non-reachable node" exception=(e, catch_backtrace())
        end
    end
    return info
end

function run_experiments(graph::PA.Graph,
                         target_nodes::Vector{Int},
                         max_depths::Vector{Int},
                         rhos::Vector{Float64},
                         gammas::Vector{Float64},
                         alphas::Vector{Float64},
                         pruning_functions::Vector{Function},
                         info_update::Function,
                         pulse_score::Function,
                         n::Int,
                         initial_bound::Bool=true)
    info = Dict(
        "elapsed_time" => Float64[],
        "preprocessing_time" => Float64[],
        "time_budget" => Float64[], 
        "alpha" => Float64[],
        "gamma" => Float64[],
        "rho" => Float64[],
        "max_depth" => Int[],
        "distance" => Int[],
    )
    for i in ProgressBar(1:length(rhos))
        for max_depth in max_depths
            for gamma in gammas
                for alpha in alphas
                    for target_node in target_nodes
                        constants = Dict("T_max" => 10.0, "alpha" => alpha, "gamma" => gamma, "rho" => rhos[i], "max_depth" => max_depth)
                        deterministic_weights = ["cost"]
                        random_weights = Dict("time" => ["mean", "variance", "covariance"])
                        prep_deterministic_weights = ["cost"]
                        prep_random_weights = Dict("time" => ["mean", "variance"])
                        parameters = Parameters(100, 
                                            false, 
                                            exploration_order, 
                                            deterministic_weights, 
                                            random_weights, 
                                            prep_deterministic_weights, 
                                            prep_random_weights)
                        graph.covariance["time"] = get_covariance_dict(graph, "time", rhos[i], max_depth)["time"]
                        instance_info = n_experiments_PaSarp(graph, target_node, parameters, constants, pruning_functions, info_update, pulse_score, initial_bound, n) 
                        append!(info["elapsed_time"], instance_info["elapsed_time"])
                        append!(info["preprocessing_time"], instance_info["preprocessing_time"])
                        append!(info["time_budget"], instance_info["time_budget"])
                        append!(info["alpha"], instance_info["alpha"])
                        append!(info["gamma"], instance_info["gamma"])
                        append!(info["rho"], instance_info["rho"])
                        append!(info["max_depth"], instance_info["max_depth"])
                        append!(info["distance"], instance_info["distance"])
                    end
                end
            end
        end
    end
    return info
end


#### UTIL ####
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