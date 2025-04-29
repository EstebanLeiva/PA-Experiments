using DataStructures: DefaultDict, PriorityQueue, enqueue!, dequeue!
using Distributions
using Random
using ProgressBars
using DataFrames
using CSV
using PulseAlgorithm
using PulseAlgorithm: Graph, Parameters, Problem, Pulse, preprocess!, run_pulse!, DefaultDict, PriorityQueue, enqueue!, dequeue!, PathInformation
const PA = PulseAlgorithm
include("util.jl")
include("sydney_loader.jl")

### Initialization ###
Random.seed!(1234)
path = raw"C:\Users\esteb\OneDrive - Universidad de los andes\Documentos\ANDES\Noveno Semestre\PA-RSPPs\PA-Experiments"

ρ = 1.0
CV = 0.8
α = 0.9
γ = 0.4
max_depth = 2
n = 1  #number of start_nodes per target_node

### Load networks ###

# Chicago Sketch (CS)
net_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_net.tntp")
flow_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_flow.tntp")

toll_factor_CS = 0.02
distance_factor_CS = 0.04

graph_CS = PA.load_graph_from_ta(net_CS, flow_CS, "CS", CV, toll_factor_CS, distance_factor_CS)
cov_CS = get_covariance_dict(graph_CS, "time", ρ, max_depth)
graph_CS.covariance["time"] = cov_CS["time"]

# Chicago Regional (CR)
net_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_net.tntp")
flow_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_flow.tntp")

toll_factor_CR = 0.1
distance_factor_CR = 0.25

graph_CR = PA.load_graph_from_ta(net_CR, flow_CR, "CR", CV, toll_factor_CR, distance_factor_CR)
cov_CR = get_covariance_dict(graph_CR, "time", ρ, max_depth)
graph_CR.covariance["time"] = cov_CR["time"]

# Sydney (SY)
net_SY = joinpath(path, raw"data\Sydney\sydney_net.tntp")
flow_SY = joinpath(path, raw"data\Sydney\sydney_flow.tntp")

toll_factor_SY = 0.1 #Use Chicago Regional's toll factor as a reference
distance_factor_SY = 0.25 #Use Chicago Regional's toll factor as a reference

graph_SY = load_graph_from_ta_without_flow_Sydney(net_SY, "SY", CV, toll_factor_SY, distance_factor_SY)
cov_SY = get_covariance_dict(graph_SY, "time", ρ, max_depth)
graph_SY.covariance["time"] = cov_SY["time"]

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
        pulse_alg.instance_info["total_length_pruned_by_feasibility"] += length(current_path)
    elseif mean > pulse_alg.problem.constants["T_max"] && prob < pulse_alg.problem.constants["alpha"]
        pass = true
        pulse_alg.instance_info["pruned_by_feasibility"] += 1
        pulse_alg.instance_info["total_length_pruned_by_feasibility"] += length(current_path)
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
        pulse_alg.instance_info["total_length_pruned_by_bounds"] += length(current_path)
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

# Vanilla
function prune_bounds_vanilla(pulse_alg::PA.Pulse, 
    current_node::Int, 
    current_path::Vector{Int},
    deterministic_info::Dict{String, Float64},
    random_info::Dict{String, Dict{String, Float64}})
   
   mean, variance, covariance_term = get_path_distribution(pulse_alg.problem.graph, current_path, pulse_alg.problem.graph.covariance["time"])
   dist = Normal(mean, √(variance + covariance_term))
   prob = cdf(dist, pulse_alg.problem.constants["T_max"])
   pass = true
   if deterministic_info["cost"] + pulse_alg.preprocessing.deterministic["cost"][current_node] <= pulse_alg.current_optimal_objective
       if current_node == pulse_alg.problem.target_node && prob >= pulse_alg.problem.constants["alpha"]
           mean, variance, covariance_term = get_path_distribution(pulse_alg.problem.graph, current_path, pulse_alg.problem.graph.covariance["time"])
           pulse_alg.current_optimal_objective = deterministic_info["cost"]
           new_path = copy(current_path)
           push!(new_path, current_node)
           pulse_alg.current_optimal_path = new_path
       end
       pass = false
       pulse_alg.instance_info["pruned_by_bounds"] += 1
       pulse_alg.instance_info["total_length_pruned_by_bounds"] += length(current_path)
   end
   #Extra prune when a path becomes infeasible
   if prob < pulse_alg.problem.constants["alpha"]
       pass = false
   end
   return pass
end

# Parameters
pruning_functions_vanilla = Vector{Function}(undef, 1)
pruning_functions_vanilla[1] = prune_bounds_vanilla
pruning_functions_pasarp = [prune_bounds, prune_feasibility]
constants = Dict("T_max" => 10.0, "alpha" => α, "gamma" => γ)
deterministic_weights = ["cost"]
random_weights = Dict("time" => ["mean", "variance", "covariance"])
prep_deterministic_weights = ["cost"]
prep_random_weights = Dict("time" => ["mean", "variance"])
params = Parameters(100, 
                    false, 
                    exploration_order, 
                    deterministic_weights, 
                    random_weights, 
                    prep_deterministic_weights, 
                    prep_random_weights)

### Computational Time Experiments ###

## CS
sampled_keys_CS = sample(collect(keys(graph_CS.nodes)), 10, replace=false)
println("Sampled keys for Chicago Sketch: ", sampled_keys_CS)

info_CS = aggregate_experiments(sampled_keys_CS, graph_CS, params, constants, pruning_functions_vanilla, info_update, pulse_score, n)
CS_nondominated_paths = info_CS["number_nondominanted_paths"]
CS_elapsed_time = info_CS["total_elapsed_time"]
CS_pruned_by_bounds = info_CS["pruned_by_bounds"]
CS_pruned_by_feasibility = info_CS["pruned_by_feasibility"]
CS_length_pruned_by_bounds = (info_CS["total_length_pruned_by_bounds"] ./ info_CS["pruned_by_bounds"])
CS_length_pruned_by_feasibility = (info_CS["total_length_pruned_by_feasibility"] ./ info_CS["pruned_by_feasibility"])
CS_preprocessing_time = info_CS["total_preprocessing_time"]

println("---------------------------------")
println("Chicago Sketch")
println("---------------------------------")

# CR
sampled_keys_CR = sample(collect(keys(graph_CR.nodes)), 10, replace=false)
println("Sampled keys for Chicago Regional: ", sampled_keys_CR)
info_CR = aggregate_experiments(sampled_keys_CR, graph_CR, params, constants, pruning_functions_vanilla, info_update, pulse_score, n)

CR_nondominated_paths = info_CR["number_nondominanted_paths"]
CR_elapsed_time = info_CR["total_elapsed_time"]
CR_pruned_by_bounds = info_CR["pruned_by_bounds"]
CR_pruned_by_feasibility = info_CR["pruned_by_feasibility"]
CR_length_pruned_by_bounds = (info_CR["total_length_pruned_by_bounds"] ./ info_CR["pruned_by_bounds"])
CR_length_pruned_by_feasibility = (info_CR["total_length_pruned_by_feasibility"] ./ info_CR["pruned_by_feasibility"])
CR_preprocessing_time = info_CR["total_preprocessing_time"]

println("---------------------------------")
println("Chicago Regional")
println("---------------------------------")


# SY
sampled_keys_SY = sample(collect(keys(graph_SY.nodes)), 10, replace=false)
println("Sampled keys for Sydney: ", sampled_keys_SY)
info_SY = aggregate_experiments(sampled_keys_SY, graph_SY, params, constants, pruning_functions_vanilla, info_update, pulse_score, n)
SY_nondominated_paths = info_SY["number_nondominanted_paths"]
SY_elapsed_time = info_SY["total_elapsed_time"]
SY_pruned_by_bounds = info_SY["pruned_by_bounds"]
SY_pruned_by_feasibility = info_SY["pruned_by_feasibility"]
SY_length_pruned_by_bounds = (info_SY["total_length_pruned_by_bounds"] ./ info_SY["pruned_by_bounds"])
SY_length_pruned_by_feasibility = (info_SY["total_length_pruned_by_feasibility"] ./ info_SY["pruned_by_feasibility"])
SY_preprocessing_time = info_SY["total_preprocessing_time"]
println("---------------------------------")
println("Sydney")
println("---------------------------------")

### Write CSV ###

df_time = DataFrame(CS_t = CS_elapsed_time, 
                    CS_t_prep = CS_preprocessing_time, 
                    CR_t = CR_elapsed_time,
                    CR_t_prep = CR_preprocessing_time,
                    SY_b = SY_elapsed_time,
                    SY_t_prep = SY_preprocessing_time)
CSV.write("time_vanilla.csv", df_time)

df_prune = DataFrame(CS_b = CS_pruned_by_bounds, 
                        CS_i = CS_pruned_by_feasibility,
                        CR_b = CR_pruned_by_bounds, 
                        CR_i = CR_pruned_by_feasibility,
                        SY_b = SY_pruned_by_bounds,
                        SY_i = SY_pruned_by_feasibility)

CSV.write("prune_vanilla.csv", df_prune)

df_length_prune = DataFrame(CS_b = CS_length_pruned_by_bounds, 
                        CS_i = CS_length_pruned_by_feasibility,
                        CR_b = CR_length_pruned_by_bounds, 
                        CR_i = CR_length_pruned_by_feasibility,
                        SY_b = SY_length_pruned_by_bounds,
                        SY_i = SY_length_pruned_by_feasibility)

CSV.write("length_prune_vanilla.csv", df_length_prune)