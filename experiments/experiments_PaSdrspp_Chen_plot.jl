using PulseAlgorithm
const PA = PulseAlgorithm
using Random
using Distributions
using DataFrames
using CSV
using DataStructures: DefaultDict, PriorityQueue, enqueue!, dequeue!
include("util_sdrspp.jl")
include("sydney_loader.jl")
include("chen_loader.jl")

Random.seed!(1234)
path = raw"C:\Users\esteb\OneDrive - Universidad de los andes\Documentos\ANDES\Noveno Semestre\PA-RSPPs\PA-Experiments"

α = 0.9
max_depth = 1
n = 10 #number of start_nodes per target_node
dominance_limit = 5

### PULSE ###
function exploration_order(pulse_alg::PA.Pulse, node::Int)
    return pulse_alg.preprocessing.random["time"]["mean"][node]
end

function pulse_score(pulse_alg::PA.Pulse, 
                     current_path::Vector{Int},
                     deterministic_info::Dict{String, Float64},
                     random_info::Dict{String, Dict{String, Float64}})
    return pulse_alg.preprocessing.random["time"]["mean"][current_path[end]]
end

function info_update(graph::PA.Graph,
                     current_node::Int, 
                     reachable_node::Int, 
                     path::Vector{Int},
                     deterministic_info::Dict{String, Float64}, 
                     random_info::Dict{String, Dict{String, Float64}})
    deterministic_info = copy(deterministic_info) #TODO: check if this copying can be done automatically in some way
    random_info = copy(random_info)

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

function prune_bounds(pulse_alg::PA.Pulse, 
    current_node::Int, 
    current_path::Vector{Int},
    deterministic_info::Dict{String, Float64},
    random_info::Dict{String, Dict{String, Float64}})
    mean = random_info["time"]["mean"] + pulse_alg.preprocessing.random["time"]["mean"][current_node]
    variance = random_info["time"]["variance"] + pulse_alg.preprocessing.random["time"]["variance"][current_node]
    dist = Normal(mean, √variance)
    prob = cdf(dist, pulse_alg.current_optimal_objective)
    if mean <= pulse_alg.current_optimal_objective && prob < pulse_alg.problem.constants["alpha"]
        pulse_alg.instance_info["pruned_by_bounds"] += 1
        pulse_alg.instance_info["total_length_pruned_by_bounds"] += length(current_path)
        return true
    end
    if mean > pulse_alg.current_optimal_objective
        pulse_alg.instance_info["pruned_by_bounds"] += 1
        pulse_alg.instance_info["total_length_pruned_by_bounds"] += length(current_path)
        return true
    end
    quant = quantile(dist, pulse_alg.problem.constants["alpha"])
    if current_node == pulse_alg.problem.target_node && quant <= pulse_alg.current_optimal_objective
        pulse_alg.current_optimal_objective = quant
        new_path = copy(current_path)
        push!(new_path, current_node)
        pulse_alg.current_optimal_path = new_path
    end
    return false
end

function prune_dominance(pulse_alg::PA.Pulse, 
       current_node::Int, 
       current_path::Vector{Int},
       deterministic_info::Dict{String, Float64},
       random_info::Dict{String, Dict{String, Float64}})
    if length(current_path) == 0
        return false
    end
    mean = random_info["time"]["mean"]
    variance = random_info["time"]["variance"] + random_info["time"]["covariance"]
    path_info = PA.PathInformation(current_path, deterministic_info, random_info)
    if (current_path[end], current_node) ∉ keys(pulse_alg.dominance)
        pulse_alg.dominance[(current_path[end], current_node)] = PriorityQueue{Tuple{Dict{String, Float64}, Dict{String, Dict{String, Float64}}}, Float64}()
        enqueue!(pulse_alg.dominance[(current_path[end], current_node)], path_info, random_info["time"]["mean"])
        return false
    end
    queue = pulse_alg.dominance[(current_path[end], current_node)]
    iter_queue = copy(queue)
    while !isempty(iter_queue)
        q_path_info = dequeue!(iter_queue)
        if q_path_info.random["time"]["mean"] < mean && q_path_info.random["time"]["variance"] < variance
            pulse_alg.instance_info["pruned_by_bounds"] += 1
            pulse_alg.instance_info["total_length_pruned_by_bounds"] += length(current_path)
            return true
        end
    end
    if length(queue) < dominance_limit
        enqueue!(queue, path_info, random_info["time"]["mean"])
    end
    return false
end

### EXPERIMENTS ###

pruning_functions = Vector{Function}()
push!(pruning_functions, prune_bounds)
pruning_functions_dominance = [prune_bounds, prune_dominance]

## Load networks ##

# Chicago Sketch (CS)
net_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_net.tntp")
flow_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_flow.tntp")
folder_CS = joinpath(path, raw"data\ChicagoSketch")

divisor_CS = 5280.0
max_speed_CS = 100.0

graph_CS = load_graph_from_ta_chen(net_CS, "")
cov_CS = get_covariance_dict_chen(graph_CS, "time", max_depth)
graph_CS.covariance["time"] = cov_CS["time"]

# Philadelphia (PH)
net_PH = joinpath(path, raw"data\Philadelphia\Philadelphia_net.tntp")
folder_PH = joinpath(path, raw"data\Philadelphia")

divisor_PH = 100.0
max_speed_PH = 100.0

graph_PH = load_graph_Philadelphia_chen(net_PH, "")
cov_PH = get_covariance_dict_chen(graph_PH, "time", max_depth)
graph_PH.covariance["time"] = cov_PH["time"]


# Sydney (SY)
net_SY = joinpath(path, raw"data\Sydney\sydney_net.tntp")
folder_SY = joinpath(path, raw"data\Sydney")

divisor_SY = 1.0
max_speed_SY = 100.0

graph_SY = load_graph_Sydney_chen(net_SY, "")
cov_SY = get_covariance_dict_chen(graph_SY, "time", max_depth)
graph_SY.covariance["time"] = cov_SY["time"]

## Computational Time Experiments ##

# CS
sampled_keys_CS = rand(collect(keys(graph_CS.nodes)), 10)
println("Sampled keys for Chicago Sketch: ", sampled_keys_CS)
pulse_info_CS, pulse_dominance_info_CS, erspa_info_CS = aggregate_experiments(sampled_keys_CS, 
                                                     graph_CS, 
                                                     α,
                                                     folder_CS,
                                                     true, 
                                                     n,
                                                     max_speed_CS,
                                                     divisor_CS,
                                                     pruning_functions,
                                                     pruning_functions_dominance,
                                                     pulse_score,
                                                     info_update)

CS_nondominated_paths_pulse = pulse_info_CS["number_nondominanted_paths"] 
CS_elapsed_time_pulse = pulse_info_CS["total_elapsed_time"] 
CS_preprocessing_time_pulse = pulse_info_CS["preprocessing_time"]
CS_pruned_by_bounds_pulse = pulse_info_CS["pruned_by_bounds"] 
CS_length_pruned_by_bounds_pulse = (pulse_info_CS["total_length_pruned_by_bounds"] ./ pulse_info_CS["pruned_by_bounds"])

CS_nondominated_paths_pulsedom = pulse_dominance_info_CS["number_nondominanted_paths"]
CS_elapsed_time_pulsedom = pulse_dominance_info_CS["total_elapsed_time"]
CS_preprocessing_time_pulsedom = pulse_dominance_info_CS["preprocessing_time"]
CS_pruned_by_bounds_pulsedom = pulse_dominance_info_CS["pruned_by_bounds"]
CS_length_pruned_by_bounds_pulsedom = (pulse_dominance_info_CS["total_length_pruned_by_bounds"] ./ pulse_dominance_info_CS["pruned_by_bounds"])
CS_pruned_by_dominance_pulsedom = pulse_dominance_info_CS["pruned_by_dominance"]
CS_length_pruned_by_dominance_pulsedom = (pulse_dominance_info_CS["total_length_pruned_by_dominance"] ./ pulse_dominance_info_CS["pruned_by_dominance"])

CS_nondominated_paths_erspa = erspa_info_CS["number_nondominanted_paths"]
CS_elapsed_time_erspa = erspa_info_CS["total_elapsed_time"]
CS_preprocessing_time_erspa = erspa_info_CS["preprocessing_time"]

println("-------------------------------------------------")
println("Chicago Sketch")
println("-------------------------------------------------")


# PH
sampled_keys_PH = sample(collect(keys(graph_PH.nodes)), 10, replace=false)
println("Sampled keys for Philadelphia: ", sampled_keys_PH)
pulse_info_PH, pulse_dominance_info_PH, erspa_info_PH = aggregate_experiments(sampled_keys_PH, 
                                                    graph_PH, 
                                                    α,
                                                    folder_PH,
                                                    true, 
                                                    n,
                                                    max_speed_PH,
                                                    divisor_PH,
                                                    pruning_functions,
                                                    pruning_functions_dominance,
                                                    pulse_score,
                                                    info_update)

PH_nondominated_paths_pulse = pulse_info_PH["number_nondominanted_paths"]
PH_elapsed_time_pulse = pulse_info_PH["total_elapsed_time"]
PH_preprocessing_time_pulse = pulse_info_PH["preprocessing_time"]
PH_pruned_by_bounds_pulse = pulse_info_PH["pruned_by_bounds"]
PH_length_pruned_by_bounds_pulse = (pulse_info_PH["total_length_pruned_by_bounds"] ./ pulse_info_PH["pruned_by_bounds"])

PH_nondominated_paths_pulsedom = pulse_dominance_info_PH["number_nondominanted_paths"]
PH_elapsed_time_pulsedom = pulse_dominance_info_PH["total_elapsed_time"]
PH_preprocessing_time_pulsedom = pulse_dominance_info_PH["preprocessing_time"]
PH_pruned_by_bounds_pulsedom = pulse_dominance_info_PH["pruned_by_bounds"]
PH_length_pruned_by_bounds_pulsedom = (pulse_dominance_info_PH["total_length_pruned_by_bounds"] ./ pulse_dominance_info_PH["pruned_by_bounds"])
PH_pruned_by_dominance_pulsedom = pulse_dominance_info_PH["pruned_by_dominance"]
PH_length_pruned_by_dominance_pulsedom = (pulse_dominance_info_PH["total_length_pruned_by_dominance"] ./ pulse_dominance_info_PH["pruned_by_dominance"])

PH_nondominated_paths_erspa = erspa_info_PH["number_nondominanted_paths"]
PH_elapsed_time_erspa = erspa_info_PH["total_elapsed_time"]
PH_preprocessing_time_erspa = erspa_info_PH["preprocessing_time"]

println("-------------------------------------------------")
println("Philadelphia")
println("-------------------------------------------------")


# SY
sampled_keys_SY = sample(collect(keys(graph_SY.nodes)), 10, replace=false)
println("Sampled keys for Sydney: ", sampled_keys_SY)
pulse_info_SY, pulse_dominance_info_SY, erspa_info_SY = aggregate_experiments(sampled_keys_SY, 
                                                    graph_SY, 
                                                    α,
                                                    folder_SY,
                                                    true, 
                                                    n,
                                                    max_speed_SY,
                                                    divisor_SY,
                                                    pruning_functions,
                                                    pruning_functions_dominance,
                                                    pulse_score,
                                                    info_update)

SY_nondominated_paths_pulse = pulse_info_SY["number_nondominanted_paths"]
SY_elapsed_time_pulse = pulse_info_SY["total_elapsed_time"]
SY_preprocessing_time_pulse = pulse_info_SY["preprocessing_time"]
SY_pruned_by_bounds_pulse = pulse_info_SY["pruned_by_bounds"]
SY_length_pruned_by_bounds_pulse = (pulse_info_SY["total_length_pruned_by_bounds"] ./ pulse_info_SY["pruned_by_bounds"])

SY_nondominated_paths_pulsedom = pulse_dominance_info_SY["number_nondominanted_paths"]
SY_elapsed_time_pulsedom = pulse_dominance_info_SY["total_elapsed_time"]
SY_preprocessing_time_pulsedom = pulse_dominance_info_SY["preprocessing_time"]
SY_pruned_by_bounds_pulsedom = pulse_dominance_info_SY["pruned_by_bounds"]
SY_length_pruned_by_bounds_pulsedom = (pulse_dominance_info_SY["total_length_pruned_by_bounds"] ./ pulse_dominance_info_SY["pruned_by_bounds"])
SY_pruned_by_dominance_pulsedom = pulse_dominance_info_SY["pruned_by_dominance"]
SY_length_pruned_by_dominance_pulsedom = (pulse_dominance_info_SY["total_length_pruned_by_dominance"] ./ pulse_dominance_info_SY["pruned_by_dominance"])

SY_nondominated_paths_erspa = erspa_info_SY["number_nondominanted_paths"]
SY_elapsed_time_erspa = erspa_info_SY["total_elapsed_time"]
SY_preprocessing_time_erspa = erspa_info_SY["preprocessing_time"]

println("-------------------------------------------------")
println("Sydney")
println("-------------------------------------------------")

### DFs ###
df_time = DataFrame(CS_pt = CS_preprocessing_time_pulse, CS_pt_dominance = CS_preprocessing_time_pulsedom, CS_pt_ERSPA = CS_preprocessing_time_erspa,
                    CS_t = CS_elapsed_time_pulse, CS_pt_dominance = CS_elapsed_time_pulsedom, CS_t_ERSPA = CS_elapsed_time_erspa, 
                    PH_pt = PH_preprocessing_time_pulse, PH_pt_dominance = PH_preprocessing_time_pulsedom, PH_pt_ERSPA = PH_preprocessing_time_erspa,
                    PH_t = PH_elapsed_time_pulse, PH_pt_dominance = PH_elapsed_time_pulsedom, PH_t_ERSPA = PH_elapsed_time_erspa,
                    SY_pt = SY_preprocessing_time_pulse, SY_pt_dominance = SY_preprocessing_time_pulsedom, SY_pt_ERSPA = SY_preprocessing_time_erspa,
                    SY_t = SY_elapsed_time_pulse, SY_pt_dominance = SY_elapsed_time_pulsedom, SY_t_ERSPA = SY_elapsed_time_erspa)

CSV.write("time_sdrspp.csv", df_time)

df_nondominated = DataFrame(CS_n = CS_nondominated_paths_pulse, CS_n_dominance = CS_nondominated_paths_pulsedom, CS_n_ERSPA = CS_nondominated_paths_erspa,
                            PH_n = PH_nondominated_paths_pulse, PH_n_dominance = PH_nondominated_paths_pulsedom, PH_n_ERSPA = PH_nondominated_paths_erspa,
                            SY_n = SY_nondominated_paths_pulse, SY_n_dominance = SY_nondominated_paths_pulsedom, SY_n_ERSPA = SY_nondominated_paths_erspa)

CSV.write("nondominated_sdrspp.csv", df_nondominated)