using ProgressBars

function get_initial_paths(graph::Graph, start_node::Int, target_node::Int, α::Float64, cov_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64})
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

function erspa_preprocess_experiments(graph::Graph, source_node::Int, target_node::Int, α::Float64, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, folder_path::String, max_speed::Float64, distance_divisor::Float64)
    file_path = joinpath(folder_path, "node_coordinates.tntp")
    node_coordinates = Vector{Tuple{Float64,Float64}}()
    open(file_path, "r") do file
        readline(file)
        for line in eachline(file)
            if occursin("Philadelphia", folder_path)
                parts = split(line, ' ')
            else
                parts = split(line, '\t')
            end
            x = parse(Float64, parts[2])
            y = parse(Float64, parts[3])
            push!(node_coordinates, (x, y))
        end
    end
    source_node = graph.nodes[source_node].name
    target_node = graph.nodes[target_node].name
    erspa = PA.initialize_ErspaStar(graph, α, covariance_dict, node_coordinates, source_node, target_node, max_speed, distance_divisor)
    PA.preprocess!(erspa)
    return erspa
end

function experiment_sdrspp(graph::Graph, source_node::Int, target_node::Int, pulse::PaSdrspp, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, α::Float64, folder_path::String, initial_bound::Bool, max_speed::Float64, distance_divisor::Float64)
    shortest_mean_path, _, _, _ = get_initial_paths(graph, source_node, target_node, α, covariance_dict)
    mean_m, variance_m, covariance_term_m = PA.get_path_distribution(graph, shortest_mean_path, covariance_dict)
    quantile_m = quantile(Normal(mean_m, √(variance_m + covariance_term_m)), α)
    if initial_bound
        elapsed_time_pulse = @elapsed begin
            PA.run_pulse(pulse, shortest_mean_path, quantile_m)
        end
    end
    erspa = erspa_preprocess_experiments(graph, source_node, target_node, α, covariance_dict, folder_path, max_speed, distance_divisor)
    elapsed_time_erspa = @elapsed begin
        PA.run_erspa(erspa)
    end
    return elapsed_time_pulse, pulse.instance_info, pulse.optimal_path, elapsed_time_erspa, erspa.instance_info, erspa.optimal_path
end

function n_experiments_sdrspp(graph::Graph, target_node::Int, α::Float64, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, folder_path::String, initial_bound::Bool, n::Int, max_speed::Float64, distance_divisor::Float64)
    pulse = PA.initialize_PaSdrspp(graph, α, covariance_dict, graph.nodes[target_node].name, graph.nodes[target_node].name, 0) ##### MODIFICATON
    PA.preprocess!(pulse)
    info_pulse = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}()
    )
    info_erspa = Dict(
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}()
    )
    source_nodes = sample(collect(keys(graph.nodes)), n, replace=false)
    for source_node in source_nodes
        shortest_mean_path, _, _, _ = get_initial_paths(graph, source_node, target_node, α, covariance_dict) ##### MODIFICATON
        pulse.max_pulse_depth = trunc(Int, length(shortest_mean_path)) ##### MODIFICATON
        pulse.source_node = source_node
        elapsed_time_pulse, instance_info_pulse, optimal_path_pulse, elapsed_time_erspa, instance_info_erspa, optimal_path_erspa = experiment_sdrspp(graph, source_node, target_node, pulse, covariance_dict, α, folder_path, initial_bound, max_speed, distance_divisor)
        if optimal_path_erspa != optimal_path_pulse
            println("-------------------------------------------------------------------")
            println("Different optimal paths")
            println("ERSPA")
            println(optimal_path_erspa)
            mean, variance, covariance = PA.get_path_distribution(graph, optimal_path_erspa, covariance_dict)
            println((mean, variance, covariance))
            println(quantile(Normal(mean, √(variance + covariance)), α))
            println("PULSE")
            println(optimal_path_pulse)
            mean, variance, covariance = PA.get_path_distribution(graph, optimal_path_pulse, covariance_dict)
            println((mean, variance, covariance))
            println(quantile(Normal(mean, √(variance + covariance)), α))
            println("-------------------------------------------------------------------")
            error("Different optimal paths")
        end
        append!(info_pulse["pruned_by_bounds"], instance_info_pulse["pruned_by_bounds"])
        append!(info_pulse["total_length_pruned_by_bounds"], instance_info_pulse["total_length_pruned_by_bounds"])
        append!(info_pulse["number_nondominanted_paths"], instance_info_pulse["number_nondominanted_paths"])
        append!(info_pulse["total_elapsed_time"], elapsed_time_pulse)

        append!(info_erspa["number_nondominanted_paths"], instance_info_erspa["number_nondominanted_paths"])
        append!(info_erspa["total_elapsed_time"], elapsed_time_erspa)

    end
    return info_pulse, info_erspa
end

function aggregate_experiments(sampled_keys::Vector{Int}, graph::Graph, α::Float64, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, folder_path::String, initial_bound::Bool, n::Int, max_speed::Float64, distance_divisor::Float64)
    total_instance_info_pulse = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}()
    )

    total_instance_info_erspa = Dict(
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}()
    )

    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info_pulse, instance_info_erspa = n_experiments_sdrspp(graph, key, α, covariance_dict, folder_path, initial_bound, n, max_speed, distance_divisor)
        
        total_instance_info_pulse["pruned_by_bounds"] = [total_instance_info_pulse["pruned_by_bounds"]; instance_info_pulse["pruned_by_bounds"]]
        total_instance_info_pulse["total_length_pruned_by_bounds"] = [total_instance_info_pulse["total_length_pruned_by_bounds"]; instance_info_pulse["total_length_pruned_by_bounds"]]
        total_instance_info_pulse["number_nondominanted_paths"] = [total_instance_info_pulse["number_nondominanted_paths"]; instance_info_pulse["number_nondominanted_paths"]]
        total_instance_info_pulse["total_elapsed_time"] = [total_instance_info_pulse["total_elapsed_time"]; instance_info_pulse["total_elapsed_time"]]

        total_instance_info_erspa["number_nondominanted_paths"] = [total_instance_info_erspa["number_nondominanted_paths"]; instance_info_erspa["number_nondominanted_paths"]]
        total_instance_info_erspa["total_elapsed_time"] = [total_instance_info_erspa["total_elapsed_time"]; instance_info_erspa["total_elapsed_time"]]
    
    end
    return total_instance_info_pulse, total_instance_info_erspa
end

### PA alone ###

function pa_experiment_sdrspp(graph::Graph, source_node::Int, target_node::Int, pulse::PaSdrspp, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, α::Float64, initial_bound::Bool)
    shortest_mean_path, _, _, _ = get_initial_paths(graph, source_node, target_node, α, covariance_dict)
    mean_m, variance_m, covariance_term_m = PA.get_path_distribution(graph, shortest_mean_path, covariance_dict)
    quantile_m = quantile(Normal(mean_m, √(variance_m + covariance_term_m)), α)
    if initial_bound
        elapsed_time_pulse = @elapsed begin
            PA.run_pulse(pulse, shortest_mean_path, quantile_m)
        end
    end
    return elapsed_time_pulse, pulse.instance_info
end

function pa_n_experiments_sdrspp(graph::Graph, target_node::Int, α::Float64, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, initial_bound::Bool, n::Int)
    pulse = PA.initialize_PaSdrspp(graph, α, covariance_dict, graph.nodes[target_node].name, graph.nodes[target_node].name)
    PA.preprocess!(pulse)
    info_pulse = Dict(
        "pruned_by_bounds" => 0,
        "total_length_pruned_by_bounds" => 0,
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0
    )
    source_nodes = sample(collect(keys(graph.nodes)), n, replace=false)
    for source_node in source_nodes
        pulse.source_node = source_node
        elapsed_time_pulse, instance_info_pulse = pa_experiment_sdrspp(graph, source_node, target_node, pulse, covariance_dict, α, initial_bound)
        info_pulse["pruned_by_bounds"] += instance_info_pulse["pruned_by_bounds"]
        info_pulse["total_length_pruned_by_bounds"] += instance_info_pulse["total_length_pruned_by_bounds"]
        info_pulse["number_nondominanted_paths"] += instance_info_pulse["number_nondominanted_paths"]
        info_pulse["total_elapsed_time"] += elapsed_time_pulse
    end
    return info_pulse
end

function pa_aggregate_experiments(sampled_keys::Vector{Int}, graph::Graph, α::Float64, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, initial_bound::Bool, n::Int)
    total_instance_info_pulse = Dict(
        "pruned_by_bounds" => 0,
        "total_length_pruned_by_bounds" => 0,
        "number_nondominanted_paths" => 0,
        "total_elapsed_time" => 0.0
    )

    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info_pulse = pa_n_experiments_sdrspp(graph, key, α, covariance_dict, initial_bound, n)
        total_instance_info_pulse["pruned_by_bounds"] += instance_info_pulse["pruned_by_bounds"]
        total_instance_info_pulse["total_length_pruned_by_bounds"] += instance_info_pulse["total_length_pruned_by_bounds"]
        total_instance_info_pulse["number_nondominanted_paths"] += instance_info_pulse["number_nondominanted_paths"]
        total_instance_info_pulse["total_elapsed_time"] += instance_info_pulse["total_elapsed_time"]
    end
    return total_instance_info_pulse
end