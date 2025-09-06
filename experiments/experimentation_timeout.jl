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

Random.seed!(1234)
path = raw"C:\Users\esteb\OneDrive - Universidad de los andes\Documentos\ANDES\Noveno Semestre\PA-RSPPs\PA-Experiments"

ρ = 1.0
CV = 0.8
α = 0.9
γ = 0.4
max_depth = 2

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

#### TIMEOUT EXPERIMENTS ####
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

## Experiment configuration ##
pruning_functions_bounds = Vector{Function}(undef, 1)
pruning_functions_bounds[1] = prune_bounds

budget = 10 * 60 #10 minutes
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
                    prep_random_weights,
                    budget,
                    false)
prep_det_weights = ["cost", "weight"]
params_weight = Parameters(100, 
                    false, 
                    exploration_order, 
                    deterministic_weights, 
                    random_weights, 
                    prep_det_weights,
                    prep_random_weights,
                    budget,
                    false)



# Chicago Regional
df = CSV.read("CR_bounds_experiment.csv", DataFrame)
timeout_df = df[df.timeout .== true, :]
println("Number of timeouts: ", nrow(timeout_df))

pair_dict_CR = Dict{Int64, Vector{Int64}}()
for row in eachrow(timeout_df)
    source = row.source
    target = row.target
    if haskey(pair_dict_CR, target)
        push!(pair_dict_CR[target], source)
    else
        pair_dict_CR[target] = [source]
    end
end

experiments_sarp(pair_dict_CR, 
                    graph_CR, 
                    params,
                    params_weight,
                    constants, 
                    pruning_functions_bounds, 
                    info_update, 
                    pulse_score, 
                    "CR_full")

# Sydney
budget = 10 * 60 #10 minutes
df = CSV.read("SY_bounds_experiment.csv", DataFrame)
timeout_df = df[df.timeout .== true, :]
println("Number of timeouts: ", nrow(timeout_df))

pair_dict_SY = Dict{Int64, Vector{Int64}}()
for row in eachrow(timeout_df)
    source = row.source
    target = row.target
    if haskey(pair_dict_SY, target)
        push!(pair_dict_SY[target], source)
    else
        pair_dict_SY[target] = [source]
    end
end

experiments_sarp(pair_dict_SY, 
                    graph_SY,
                    params,
                    params_weight,
                    constants, 
                    pruning_functions_bounds, 
                    info_update, 
                    pulse_score, 
                    "SY_timeout")