using ProgressBars
include("erspa_star.jl")

### UTIL ###

function get_initial_paths(graph::PA.Graph, start_node::Int, target_node::Int, α::Float64, cov_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64})
    shortest_mean_path = PA.dijkstra(graph, start_node, target_node, "time", "mean")
    shortest_var_path = PA.dijkstra(graph, start_node, target_node, "time", "variance")
    return shortest_mean_path, shortest_var_path
end

function erspa_preprocess_experiments(graph::PA.Graph, source_node::Int, target_node::Int, α::Float64, covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, folder_path::String, max_speed::Float64, distance_divisor::Float64)
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
    erspa = initialize_ErspaStar(graph, α, covariance_dict, node_coordinates, source_node, target_node, max_speed, distance_divisor)
    preprocessing_time = @elapsed begin
        preprocess_erspa!(erspa)
    end
    return erspa, preprocessing_time
end

function experiment_sdrspp(graph::PA.Graph, source_node::Int, target_node::Int, 
                           covariance_dict::PA.DefaultDict{Tuple{Int,Int,Int,Int},Float64}, α::Float64, 
                           folder_path::String, initial_bound::Bool, max_speed::Float64, distance_divisor::Float64, 
                           parameters::PA.Parameters, problem::PA.Problem, pruning_functions::Vector{Function},
                           pruning_functions_dominance::Vector{Function}, pulse_score::Function, info_update::Function)
    experiment_info = Dict(
        "instance_info_pulse" => Dict(
            "pruned_by_bounds" => 0,
            "total_length_pruned_by_bounds" => 0,
            "number_nondominanted_paths" => 0,
            "total_elapsed_time" => 0.0, 
            "preprocessing_time" => 0.0,
            "optimal_path" => Vector{Int}(),
        ),
        "instance_info_pulse_dominance" => Dict(
            "pruned_by_bounds" => 0,
            "pruned_by_dominance" => 0,
            "total_length_pruned_by_dominance" => 0,
            "total_length_pruned_by_bounds" => 0,
            "number_nondominanted_paths" => 0,
            "total_elapsed_time" => 0.0,
            "preprocessing_time" => 0.0,
            "optimal_path" => Vector{Int}(),
        ),
        "instance_info_erspa" => Dict(
            "number_nondominanted_paths" => 0,
            "total_elapsed_time" => 0.0, 
            "preprocessing_time" => 0.0,
            "optimal_path" => Vector{Int}(),
        )
    )
    
    shortest_mean_path, _ = get_initial_paths(graph, source_node, target_node, α, covariance_dict)
    mean_m, variance_m, covariance_term_m = get_path_distribution(graph, shortest_mean_path, covariance_dict)
    quantile_m = quantile(Normal(mean_m, √(variance_m + covariance_term_m)), α)

    # Pulse
    pulse = PA.Pulse(problem, parameters)
    pulse.parameters.max_pulse_depth = length(shortest_mean_path)
    pulse.instance_info["pruned_by_bounds"] = 0
    pulse.instance_info["pruned_by_feasibility"] = 0
    pulse.instance_info["total_length_pruned_by_bounds"] = 0
    pulse.instance_info["total_length_pruned_by_feasibility"] = 0
    pulse.instance_info["number_nondominanted_paths"] = 0

    preprocessing_time = @elapsed begin
        PA.preprocess!(pulse)
    end
    # Pulse normal
    elapsed_time_pulse = @elapsed begin
        PA.run_pulse!(pulse, info_update, pruning_functions, pulse_score, shortest_mean_path, quantile_m)
    end
    experiment_info["instance_info_pulse"]["optimal_path"] = pulse.optimal_path
    experiment_info["instance_info_pulse"]["pruned_by_bounds"] = pulse.instance_info["pruned_by_bounds"]
    experiment_info["instance_info_pulse"]["total_length_pruned_by_bounds"] = pulse.instance_info["total_length_pruned_by_bounds"]
    experiment_info["instance_info_pulse"]["number_nondominanted_paths"] = pulse.instance_info["number_nondominanted_paths"]
    experiment_info["instance_info_pulse"]["total_elapsed_time"] = elapsed_time_pulse
    experiment_info["instance_info_pulse"]["preprocessing_time"] = preprocessing_time

    # Pulse Dominance
    pulse = PA.Pulse(problem, parameters)
    pulse.parameters.max_pulse_depth = length(shortest_mean_path)
    pulse.instance_info["pruned_by_bounds"] = 0
    pulse.instance_info["pruned_by_dominance"] = 0
    pulse.instance_info["total_length_pruned_by_bounds"] = 0
    pulse.instance_info["total_length_pruned_by_dominance"] = 0
    pulse.instance_info["number_nondominanted_paths"] = 0
    preprocessing_time = @elapsed begin
        PA.preprocess!(pulse)
    end
    elapsed_time_pulse_dominance = @elapsed begin
        PA.run_pulse!(pulse, info_update, pruning_functions_dominance, pulse_score, shortest_mean_path, quantile_m)
    end

    experiment_info["instance_info_pulse_dominance"]["optimal_path"] = pulse.optimal_path
    experiment_info["instance_info_pulse_dominance"]["pruned_by_bounds"] = pulse.instance_info["pruned_by_bounds"]
    experiment_info["instance_info_pulse_dominance"]["pruned_by_dominance"] = pulse.instance_info["pruned_by_dominance"]
    experiment_info["instance_info_pulse_dominance"]["total_length_pruned_by_dominance"] = pulse.instance_info["total_length_pruned_by_dominance"]
    experiment_info["instance_info_pulse_dominance"]["total_length_pruned_by_bounds"] = pulse.instance_info["total_length_pruned_by_bounds"]
    experiment_info["instance_info_pulse_dominance"]["number_nondominanted_paths"] = pulse.instance_info["number_nondominanted_paths"]
    experiment_info["instance_info_pulse_dominance"]["total_elapsed_time"] = elapsed_time_pulse_dominance
    experiment_info["instance_info_pulse_dominance"]["preprocessing_time"] = preprocessing_time

    erspa, preprocessing_time_erspa = erspa_preprocess_experiments(graph, source_node, target_node, α, covariance_dict, folder_path, max_speed, distance_divisor)
    elapsed_time_erspa = @elapsed begin
        run_erspa(erspa)
    end
    experiment_info["instance_info_erspa"]["optimal_path"] = erspa.optimal_path
    experiment_info["instance_info_erspa"]["number_nondominanted_paths"] = length(erspa.SE)
    experiment_info["instance_info_erspa"]["total_elapsed_time"] = elapsed_time_erspa
    experiment_info["instance_info_erspa"]["preprocessing_time"] = preprocessing_time_erspa
    return experiment_info
end

function n_experiments_sdrspp(graph::PA.Graph, target_node::Int, α::Float64, folder_path::String, 
                              initial_bound::Bool, n::Int, max_speed::Float64, distance_divisor::Float64,
                              pruning_functions::Vector{Function}, pruning_functions_dominance::Vector{Function},
                              pulse_score::Function, info_update::Function)
    info_pulse = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(),
        "preprocessing_time" => Vector{Float64}(),
    )

    info_pulse_dominance = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "pruned_by_dominance" => Vector{Int}(),
        "total_length_pruned_by_dominance" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(),
        "preprocessing_time" => Vector{Float64}(),
    )

    info_erspa = Dict(
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(), 
        "preprocessing_time" => Vector{Float64}(),
    )
    source_nodes = sample(collect(keys(graph.nodes)), n, replace=false)
    for source_node in source_nodes
        covariance_dict = graph.covariance["time"]
        shortest_mean_path, _ = get_initial_paths(graph, source_node, target_node, α, covariance_dict)
        
        # Pulse parameters
        deterministic_weights = ["cost"]
        random_weights = Dict("time" => ["mean", "variance", "covariance"])
        prep_deterministic_weights = []
        prep_random_weights = Dict("time" => ["mean", "variance"])
        params = PA.Parameters(trunc(Int, length(shortest_mean_path)), 
                            false, 
                            exploration_order, 
                            deterministic_weights, 
                            random_weights, 
                            prep_deterministic_weights, 
                            prep_random_weights)

        constants = Dict("alpha" => α)
        problem = PA.Problem(graph, source_node, target_node, true, constants)

        try 
            experiment_info = experiment_sdrspp(graph, source_node, target_node, covariance_dict, α, folder_path, initial_bound, 
                                            max_speed, distance_divisor, params, problem, pruning_functions, 
                                            pruning_functions_dominance, pulse_score, info_update)
            instance_info_pulse = experiment_info["instance_info_pulse"]
            instance_info_pulse_dominance = experiment_info["instance_info_pulse_dominance"]
            instance_info_erspa = experiment_info["instance_info_erspa"]
    
            if instance_info_erspa["optimal_path"] != instance_info_pulse["optimal_path"] && instance_info_pulse_dominance["optimal_path"] != instance_info_pulse["optimal_path"]
                println("-------------------------------------------------------------------")
                println("Different optimal paths")
                println("-------------------------------------------------------------------")
    
                println("ERSPA")
                optimal_path_erspa = instance_info_erspa["optimal_path"]
                println(optimal_path_erspa)
                mean, variance, covariance = get_path_distribution(graph, optimal_path_erspa, covariance_dict)
                println((mean, variance, covariance))
                println(quantile(Normal(mean, √(variance + covariance)), α))
                println("-------------------------------------------------------------------")
    
                println("PULSE")
                optimal_path_pulse = instance_info_pulse["optimal_path"]
                println(optimal_path_pulse)
                mean, variance, covariance = get_path_distribution(graph, optimal_path_pulse, covariance_dict)
                println((mean, variance, covariance))
                println(quantile(Normal(mean, √(variance + covariance)), α))
                println("-------------------------------------------------------------------")
    
                println("PULSE Dominance")
                optimal_path_pulse_dominance = instance_info_pulse_dominance["optimal_path"]
                println(optimal_path_pulse_dominance)
                mean, variance, covariance = get_path_distribution(graph, optimal_path_pulse_dominance, covariance_dict)
                println((mean, variance, covariance))
                println(quantile(Normal(mean, √(variance + covariance)), α))
                println("-------------------------------------------------------------------")
                error("Different optimal paths")
            end
    
            append!(info_pulse["pruned_by_bounds"], instance_info_pulse["pruned_by_bounds"])
            append!(info_pulse["total_length_pruned_by_bounds"], instance_info_pulse["total_length_pruned_by_bounds"])
            append!(info_pulse["number_nondominanted_paths"], instance_info_pulse["number_nondominanted_paths"])
            append!(info_pulse["total_elapsed_time"], instance_info_pulse["total_elapsed_time"])
            append!(info_pulse["preprocessing_time"], instance_info_pulse["preprocessing_time"])
    
            append!(info_pulse_dominance["pruned_by_bounds"], instance_info_pulse_dominance["pruned_by_bounds"])
            append!(info_pulse_dominance["pruned_by_dominance"], instance_info_pulse_dominance["pruned_by_dominance"])
            append!(info_pulse_dominance["total_length_pruned_by_dominance"], instance_info_pulse_dominance["total_length_pruned_by_dominance"])
            append!(info_pulse_dominance["total_length_pruned_by_bounds"], instance_info_pulse_dominance["total_length_pruned_by_bounds"])
            append!(info_pulse_dominance["number_nondominanted_paths"], instance_info_pulse_dominance["number_nondominanted_paths"])
            append!(info_pulse_dominance["total_elapsed_time"], instance_info_pulse_dominance["total_elapsed_time"])
            append!(info_pulse_dominance["preprocessing_time"], instance_info_pulse_dominance["preprocessing_time"])
    
            append!(info_erspa["number_nondominanted_paths"], instance_info_erspa["number_nondominanted_paths"])
            append!(info_erspa["total_elapsed_time"], instance_info_erspa["total_elapsed_time"])
            append!(info_erspa["preprocessing_time"], instance_info_erspa["preprocessing_time"])
        catch e
            @warn "Error in experiment_sdrspp" exception=(e, catch_backtrace())
        end
    end
    return info_pulse, info_pulse_dominance, info_erspa
end

function aggregate_experiments(sampled_keys::Vector{Int}, graph::PA.Graph, α::Float64, folder_path::String, 
                                initial_bound::Bool, n::Int, max_speed::Float64, distance_divisor::Float64, 
                                pruning_functions::Vector{Function}, pruning_functions_dominance::Vector{Function},
                                pulse_score::Function, info_update::Function)
    total_instance_info_pulse = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(),
        "preprocessing_time" => Vector{Float64}()
    )

    total_instance_info_pulse_dominance = Dict(
        "pruned_by_bounds" => Vector{Int}(),
        "pruned_by_dominance" => Vector{Int}(),
        "total_length_pruned_by_dominance" => Vector{Int}(),
        "total_length_pruned_by_bounds" => Vector{Int}(),
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(),
        "preprocessing_time" => Vector{Float64}()
    )

    total_instance_info_erspa = Dict(
        "number_nondominanted_paths" => Vector{Int}(),
        "total_elapsed_time" => Vector{Float64}(), 
        "preprocessing_time" => Vector{Float64}()
    )

    for i in ProgressBar(1:length(sampled_keys))
        key = sampled_keys[i]
        instance_info_pulse, instance_info_pulse_dominance, instance_info_erspa = n_experiments_sdrspp(graph, key, α, folder_path, initial_bound, n, max_speed, distance_divisor, pruning_functions, pruning_functions_dominance, pulse_score, info_update)
        
        total_instance_info_pulse["pruned_by_bounds"] = [total_instance_info_pulse["pruned_by_bounds"]; instance_info_pulse["pruned_by_bounds"]]
        total_instance_info_pulse["total_length_pruned_by_bounds"] = [total_instance_info_pulse["total_length_pruned_by_bounds"]; instance_info_pulse["total_length_pruned_by_bounds"]]
        total_instance_info_pulse["number_nondominanted_paths"] = [total_instance_info_pulse["number_nondominanted_paths"]; instance_info_pulse["number_nondominanted_paths"]]
        total_instance_info_pulse["total_elapsed_time"] = [total_instance_info_pulse["total_elapsed_time"]; instance_info_pulse["total_elapsed_time"]]
        total_instance_info_pulse["preprocessing_time"] = [total_instance_info_pulse["preprocessing_time"]; instance_info_pulse["preprocessing_time"]]

        total_instance_info_pulse_dominance["pruned_by_bounds"] = [total_instance_info_pulse_dominance["pruned_by_bounds"]; instance_info_pulse_dominance["pruned_by_bounds"]]
        total_instance_info_pulse_dominance["pruned_by_dominance"] = [total_instance_info_pulse_dominance["pruned_by_dominance"]; instance_info_pulse_dominance["pruned_by_dominance"]]
        total_instance_info_pulse_dominance["total_length_pruned_by_dominance"] = [total_instance_info_pulse_dominance["total_length_pruned_by_dominance"]; instance_info_pulse_dominance["total_length_pruned_by_dominance"]]
        total_instance_info_pulse_dominance["total_length_pruned_by_bounds"] = [total_instance_info_pulse_dominance["total_length_pruned_by_bounds"]; instance_info_pulse_dominance["total_length_pruned_by_bounds"]]
        total_instance_info_pulse_dominance["number_nondominanted_paths"] = [total_instance_info_pulse_dominance["number_nondominanted_paths"]; instance_info_pulse_dominance["number_nondominanted_paths"]]
        total_instance_info_pulse_dominance["total_elapsed_time"] = [total_instance_info_pulse_dominance["total_elapsed_time"]; instance_info_pulse_dominance["total_elapsed_time"]]
        total_instance_info_pulse_dominance["preprocessing_time"] = [total_instance_info_pulse_dominance["preprocessing_time"]; instance_info_pulse_dominance["preprocessing_time"]]

        total_instance_info_erspa["number_nondominanted_paths"] = [total_instance_info_erspa["number_nondominanted_paths"]; instance_info_erspa["number_nondominanted_paths"]]
        total_instance_info_erspa["total_elapsed_time"] = [total_instance_info_erspa["total_elapsed_time"]; instance_info_erspa["total_elapsed_time"]]
        total_instance_info_erspa["preprocessing_time"] = [total_instance_info_erspa["preprocessing_time"]; instance_info_erspa["preprocessing_time"]]
    
    end
    return total_instance_info_pulse, total_instance_info_pulse_dominance, total_instance_info_erspa
end