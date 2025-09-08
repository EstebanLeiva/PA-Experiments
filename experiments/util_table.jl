using DataStructures: DefaultDict, PriorityQueue, enqueue!, dequeue!
using Distributions
using Random
using ProgressBars
using DataFrames
using CSV
using PulseAlgorithm
using PulseAlgorithm: Graph, Parameters, Problem, Pulse, preprocess!, run_pulse!, DefaultDict, PriorityQueue, enqueue!, dequeue!, PathInformation
const PA = PulseAlgorithm
include("util_new.jl")
include("sydney_loader.jl")

### Define Algorithms ###

# PaSarp
function info_update(graph::Graph,
    current_node::Int, 
    reachable_node::Int, 
    path::Vector{Int},
    deterministic_info::Dict{String, Float64}, 
    random_info::Dict{String, Dict{String, Float64}})

    deterministic_info["cost"] += graph.nodes[current_node].links[reachable_node].deterministic["cost"]
    random_info["time"]["mean"] += graph.nodes[current_node].links[reachable_node].random["time"]["mean"]
    random_info["time"]["variance"] += graph.nodes[current_node].links[reachable_node].random["time"]["variance"]
    n = length(path)
    if n > 1
        last_node = path[end]
        current_node = reachable_node
        covariance_sum = 0.0
        for i in 1:n-1
            covariance_sum += 2 * graph.covariance["time"][(path[i], path[i+1], last_node, current_node)]
        end
    else
        covariance_sum =  0.0
    end
    random_info["time"]["covariance"] += covariance_sum
    return deterministic_info, random_info
end

function prune_feasibility(pulse_alg::Pulse, 
          current_node::Int, 
          current_path::Vector{Int},
          deterministic_info::Dict{String, Float64},
          random_info::Dict{String, Dict{String, Float64}})
    
    pass = false
    mean = random_info["time"]["mean"] + pulse_alg.preprocessing.random["time"]["mean"][current_node]
    variance = random_info["time"]["variance"] + pulse_alg.preprocessing.random["time"]["variance"][current_node]
    dist = Normal(mean, √variance)
    prob = cdf(dist, pulse_alg.problem.constants["T_max"])
    if mean > pulse_alg.problem.constants["T_max"]
        pass = true
        pulse_alg.instance_info["pruned_by_feasibility"] += 1
        pulse_alg.instance_info["length_pruned_by_feasibility"] += length(current_path)
    elseif mean > pulse_alg.problem.constants["T_max"] && prob < pulse_alg.problem.constants["alpha"]
        pass = true
        pulse_alg.instance_info["pruned_by_feasibility"] += 1
        pulse_alg.instance_info["length_pruned_by_feasibility"] += length(current_path)
    end
    return pass
end

function prune_bounds(pulse_alg::Pulse, 
     current_node::Int, 
     current_path::Vector{Int},
     deterministic_info::Dict{String, Float64},
     random_info::Dict{String, Dict{String, Float64}})
    
    pass = true
    if deterministic_info["cost"] + pulse_alg.preprocessing.deterministic["cost"][current_node] <= pulse_alg.current_optimal_objective
        if current_node == pulse_alg.problem.target_node
            pulse_alg.current_optimal_objective = deterministic_info["cost"]
            new_path = copy(current_path)
            push!(new_path, current_node)
            pulse_alg.current_optimal_path = new_path
        end
        pass = false
        pulse_alg.instance_info["pruned_by_bounds"] += 1
        pulse_alg.instance_info["length_pruned_by_bounds"] += length(current_path)
    end
    return pass
end

function exploration_order(pulse_alg::Pulse, node::Int)
    return pulse_alg.preprocessing.deterministic["cost"][node]
end

function pulse_score(pulse_alg::Pulse, 
    current_path::Vector{Int},
    deterministic_info::Dict{String, Float64},
    random_info::Dict{String, Dict{String, Float64}})
    return pulse_alg.preprocessing.deterministic["cost"][current_path[end]]
end

function load_graphs(rho, max_depth)
    net_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_net.tntp")
    flow_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_flow.tntp")

    toll_factor_CS = 0.02
    distance_factor_CS = 0.04

    graph_CS = PA.load_graph_from_ta(net_CS, flow_CS, "CS", CV, toll_factor_CS, distance_factor_CS)
    cov_CS = get_covariance_dict(graph_CS, "time", rho, max_depth)
    graph_CS.covariance["time"] = cov_CS["time"]

    # Chicago Regional (CR)
    net_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_net.tntp")
    flow_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_flow.tntp")

    toll_factor_CR = 0.1
    distance_factor_CR = 0.25

    graph_CR = PA.load_graph_from_ta(net_CR, flow_CR, "CR", CV, toll_factor_CR, distance_factor_CR)
    cov_CR = get_covariance_dict(graph_CR, "time", rho, max_depth)
    graph_CR.covariance["time"] = cov_CR["time"]

    # Sydney (SY)
    net_SY = joinpath(path, raw"data\Sydney\sydney_net.tntp")
    flow_SY = joinpath(path, raw"data\Sydney\sydney_flow.tntp")

    toll_factor_SY = 0.1 #Use Chicago Regional's toll factor as a reference
    distance_factor_SY = 0.25 #Use Chicago Regional's toll factor as a reference

    graph_SY = load_graph_from_ta_without_flow_Sydney(net_SY, "SY", CV, toll_factor_SY, distance_factor_SY)
    cov_SY = get_covariance_dict(graph_SY, "time", rho, max_depth)
    graph_SY.covariance["time"] = cov_SY["time"]

    return (graph_CS, graph_CR, graph_SY)
end

function sarp_table(pulse, info_update, pruning_functions, pulse_score)
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

    return elapsed_time, preprocessing_time, pulse.parameters.timeout
end

function time_sarp(graph, source, target, alpha, gamma, multiplier, reference_time_budget)
    pruning_functions = [prune_feasibility, prune_bounds]
    time_budget = reference_time_budget * multiplier
    constants = Dict("T_max" => time_budget, "alpha" => alpha, "gamma" => gamma)
    deterministic_weights = ["cost"]
    random_weights = Dict("time" => ["mean", "variance", "covariance"])
    prep_deterministic_weights = ["cost"]
    prep_random_weights = Dict("time" => ["mean", "variance"])
    timeout_budget = 1 * 60 # 1 minute
    params = Parameters(100, 
                        false, 
                        exploration_order, 
                        deterministic_weights, 
                        random_weights, 
                        prep_deterministic_weights,
                        prep_random_weights,
                        timeout_budget,
                        false)

    prep_det_weights = ["cost", "weight"]
    params_weight = Parameters(100, 
                    false, 
                    exploration_order, 
                    deterministic_weights, 
                    random_weights, 
                    prep_det_weights,
                    prep_random_weights,
                    timeout_budget,
                    false)
    problem = PA.Problem(graph, source, target, true, constants)
    pulse_weight = PA.Pulse(problem, params_weight)
    pulse = PA.Pulse(problem, params)
    pulse.instance_info["pruned_by_bounds"] = 0
    pulse.instance_info["pruned_by_feasibility"] = 0
    pulse.instance_info["length_pruned_by_bounds"] = 0
    pulse.instance_info["length_pruned_by_feasibility"] = 0
    pulse.instance_info["number_nondominated_paths"] = 0
    elapsed_time = Inf
    distance = 0
    timeout = false
    try
        elapsed_time, preprocessing_time, timeout = sarp_table(pulse, info_update, pruning_functions, pulse_score)
    catch e
        @warn "Non-reachable node" exception=(e, catch_backtrace())
    end

    distance = length(PA.dijkstra(graph, source, target, "weight"))

    return (elapsed_time, time_budget, distance, timeout)

end

function experiment_table(ct_alpha, ct_rho, ct_gamma, ct_max_depth)
    results = DataFrame(
        Network = String[],
        Rho = Float64[],
        Time_Budget = Float64[],
        Max_Depth = Int[],
        Distance = Int[],
        Timeout = Bool[],
        Elapsed_Time = Float64[]
    )

    graph_CS, graph_CR, graph_SY = load_graphs(ct_rho, ct_max_depth)
    od_pairs_CS = sample_od_pairs(graph_CS, 1, 1)
    od_pairs_CR = sample_od_pairs(graph_CR, 1, 1)
    od_pairs_SY = sample_od_pairs(graph_SY, 1, 1)

    # get the first key of od_pairs_CS
    (source_CS, target) = first(od_pairs_CS)
    target_CS = target[1]
    (source_CR, target_CR) = first(od_pairs_CR)
    target_CR = target_CR[1]
    (source_SY, target_SY) = first(od_pairs_SY)
    target_SY = target_SY[1]

    reference_time_budget_CS, _ = get_timeBudget(graph_CS, source_CS, target_CS, ct_alpha, ct_gamma, graph_CS.covariance["time"])
    reference_time_budget_CR, _ = get_timeBudget(graph_CR, source_CR, target_CR, ct_alpha, ct_gamma, graph_CR.covariance["time"])
    reference_time_budget_SY, _ = get_timeBudget(graph_SY, source_SY, target_SY, ct_alpha, ct_gamma, graph_SY.covariance["time"])

    println("Reference Time Budgets:")
    println("CS: ", reference_time_budget_CS)
    println("CR: ", reference_time_budget_CR)
    println("SY: ", reference_time_budget_SY)

    ct_multiplier = 1
    budget_multipliers = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    rhos = [0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0]
    max_depths = [1, 2, 3]

    for multiplier in budget_multipliers
        graph_CS, graph_CR, graph_SY = load_graphs(ct_rho, ct_max_depth)
        for (source, target) in od_pairs_CS
            target = target[1]
            elapsed_time_CS, time_budget_CS, distance_CS, timeout_CS = time_sarp(graph_CS, source, target, ct_alpha, ct_gamma, multiplier, reference_time_budget_CS)
            push!(results, (Network = "CS", Rho = ct_rho, Time_Budget = time_budget_CS, Max_Depth = ct_max_depth, Distance = distance_CS, Timeout = timeout_CS, Elapsed_Time = elapsed_time_CS))
        end
        for (source, target) in od_pairs_CR
            target = target[1]
            elapsed_time_CR, time_budget_CR, distance_CR, timeout_CR = time_sarp(graph_CR, source, target, ct_alpha, ct_gamma, multiplier, reference_time_budget_CR)
            push!(results, (Network = "CR", Rho = ct_rho, Time_Budget = time_budget_CR, Max_Depth = ct_max_depth, Distance = distance_CR, Timeout = timeout_CR, Elapsed_Time = elapsed_time_CR))
        end
        for (source, target) in od_pairs_SY
            target = target[1]
            elapsed_time_SY, time_budget_SY, distance_SY, timeout_SY = time_sarp(graph_SY, source, target, ct_alpha, ct_gamma, multiplier, reference_time_budget_SY)
            push!(results, (Network = "SY", Rho = ct_rho, Time_Budget = time_budget_SY, Max_Depth = ct_max_depth, Distance = distance_SY, Timeout = timeout_SY, Elapsed_Time = elapsed_time_SY))
        end
    end

    for rho in rhos
       graph_CS, graph_CR, graph_SY = load_graphs(rho, ct_max_depth)
       for (source, target) in od_pairs_CS
           target = target[1]
           elapsed_time_CS, time_budget_CS, distance_CS, timeout_CS = time_sarp(graph_CS, source, target, ct_alpha, ct_gamma, ct_multiplier, reference_time_budget_CS)
           push!(results, (Network = "CS", Rho = rho, Time_Budget = time_budget_CS, Max_Depth = ct_max_depth, Distance = distance_CS, Timeout = timeout_CS, Elapsed_Time = elapsed_time_CS))
       end
       for (source, target) in od_pairs_CR
           target = target[1]
           elapsed_time_CR, time_budget_CR, distance_CR, timeout_CR = time_sarp(graph_CR, source, target, ct_alpha, ct_gamma, ct_multiplier, reference_time_budget_CR)
           push!(results, (Network = "CR", Rho = rho, Time_Budget = time_budget_CR, Max_Depth = ct_max_depth, Distance = distance_CR, Timeout = timeout_CR, Elapsed_Time = elapsed_time_CR))
       end
       for (source, target) in od_pairs_SY
           target = target[1]
           elapsed_time_SY, time_budget_SY, distance_SY, timeout_SY = time_sarp(graph_SY, source, target, ct_alpha, ct_gamma, ct_multiplier, reference_time_budget_SY)
           push!(results, (Network = "SY", Rho = rho, Time_Budget = time_budget_SY, Max_Depth = ct_max_depth, Distance = distance_SY, Timeout = timeout_SY, Elapsed_Time = elapsed_time_SY))
       end
    end

    for max_depth in max_depths
       graph_CS, graph_CR, graph_SY = load_graphs(ct_rho, max_depth)
       for (source, target) in od_pairs_CS
           target = target[1]
           elapsed_time_CS, time_budget_CS, distance_CS, timeout_CS = time_sarp(graph_CS, source, target, ct_alpha, ct_gamma, ct_multiplier, reference_time_budget_CS)
           push!(results, (Network = "CS", Rho = ct_rho, Time_Budget = time_budget_CS, Max_Depth = max_depth, Distance = distance_CS, Timeout = timeout_CS, Elapsed_Time = elapsed_time_CS))
       end
       for (source, target) in od_pairs_CR
           target = target[1]
           elapsed_time_CR, time_budget_CR, distance_CR, timeout_CR = time_sarp(graph_CR, source, target, ct_alpha, ct_gamma, ct_multiplier, reference_time_budget_CR)
           push!(results, (Network = "CR", Rho = ct_rho, Time_Budget = time_budget_CR, Max_Depth = max_depth, Distance = distance_CR, Timeout = timeout_CR, Elapsed_Time = elapsed_time_CR))
       end
       for (source, target) in od_pairs_SY
           target = target[1]
           elapsed_time_SY, time_budget_SY, distance_SY, timeout_SY = time_sarp(graph_SY, source, target, ct_alpha, ct_gamma, ct_multiplier, reference_time_budget_SY)
           push!(results, (Network = "SY", Rho = ct_rho, Time_Budget = time_budget_SY, Max_Depth = max_depth, Distance = distance_SY, Timeout = timeout_SY, Elapsed_Time = elapsed_time_SY))
       end
    end

    CSV.write("table.csv", results)

    return results
end