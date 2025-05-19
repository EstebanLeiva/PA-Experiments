using DataStructures: DefaultDict, PriorityQueue, enqueue!, dequeue!
using Distributions
using Random
using ProgressBars
using DataFrames
using CSV
using PulseAlgorithm
using PulseAlgorithm: Graph, Parameters, Problem, Pulse, preprocess!, run_pulse!, DefaultDict, PriorityQueue, enqueue!, dequeue!, PathInformation
const PA = PulseAlgorithm
include("util_regression.jl")
include("sydney_loader.jl")

### Initialization ###
Random.seed!(1234)
path = raw"C:\Users\esteb\OneDrive - Universidad de los andes\Documentos\ANDES\Noveno Semestre\PA-RSPPs\PA-Experiments"

### Load networks ###
CV = 0.8

# Chicago Sketch (CS)
net_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_net.tntp")
flow_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_flow.tntp")

toll_factor_CS = 0.02
distance_factor_CS = 0.04

graph_CS = PA.load_graph_from_ta(net_CS, flow_CS, "CS", CV, toll_factor_CS, distance_factor_CS)


# Chicago Regional (CR)
net_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_net.tntp")
flow_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_flow.tntp")

toll_factor_CR = 0.1
distance_factor_CR = 0.25

graph_CR = PA.load_graph_from_ta(net_CR, flow_CR, "CR", CV, toll_factor_CR, distance_factor_CR)

# Sydney (SY)
net_SY = joinpath(path, raw"data\Sydney\sydney_net.tntp")
flow_SY = joinpath(path, raw"data\Sydney\sydney_flow.tntp")

toll_factor_SY = 0.1 #Use Chicago Regional's toll factor as a reference
distance_factor_SY = 0.25 #Use Chicago Regional's toll factor as a reference

graph_SY = load_graph_from_ta_without_flow_Sydney(net_SY, "SY", CV, toll_factor_SY, distance_factor_SY)

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
    elseif mean > pulse_alg.problem.constants["T_max"] && prob < pulse_alg.problem.constants["alpha"]
        pass = true
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

### Regression Experiments ###
sn = 2
n = 2
alphas = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
gammas = [0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0]
rhos = [0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0]
max_depths = [1, 2, 3]

sampled_keys_CS = sample(collect(keys(graph_CS.nodes)), sn, replace=false)
info_CS = run_experiments(graph_CS, sampled_keys_CS, max_depths, rhos, gammas, alphas, [prune_feasibility, prune_bounds], info_update, pulse_score, n)

sampled_keys_CR = sample(collect(keys(graph_CR.nodes)), sn, replace=false)
info_CR = run_experiments(graph_CR, sampled_keys_CR, max_depths, rhos, gammas, alphas, [prune_feasibility, prune_bounds], info_update, pulse_score, n)

sampled_keys_SY = sample(collect(keys(graph_SY.nodes)), sn, replace=false)
info_SY = run_experiments(graph_SY, sampled_keys_SY, max_depths, rhos, gammas, alphas, [prune_feasibility, prune_bounds], info_update, pulse_score, n)

CS_elapsed_time = info_CS["elapsed_time"]
CS_preprocessing_time = info_CS["preprocessing_time"]
CS_time_budget = info_CS["time_budget"]
CS_alpha = info_CS["alpha"]
CS_gamma = info_CS["gamma"]
CS_rho = info_CS["rho"]
CS_max_depth = info_CS["max_depth"]

CR_elapsed_time = info_CR["elapsed_time"]
CR_preprocessing_time = info_CR["preprocessing_time"]
CR_time_budget = info_CR["time_budget"]
CR_alpha = info_CR["alpha"]
CR_gamma = info_CR["gamma"]
CR_rho = info_CR["rho"]
CR_max_depth = info_CR["max_depth"]

SY_elapsed_time = info_SY["elapsed_time"]
SY_preprocessing_time = info_SY["preprocessing_time"]
SY_time_budget = info_SY["time_budget"]
SY_alpha = info_SY["alpha"]
SY_gamma = info_SY["gamma"]
SY_rho = info_SY["rho"]
SY_max_depth = info_SY["max_depth"]

df_reg = DataFrame(
    Elapsed_Time = vcat(CS_elapsed_time, CR_elapsed_time, SY_elapsed_time),
    Preprocessing_Time = vcat(CS_preprocessing_time, CR_preprocessing_time, SY_preprocessing_time),
    Time_Budget = vcat(CS_time_budget, CR_time_budget, SY_time_budget),
    Alpha = vcat(CS_alpha, CR_alpha, SY_alpha),
    Gamma = vcat(CS_gamma, CR_gamma, SY_gamma),
    Rho = vcat(CS_rho, CR_rho, SY_rho),
    Max_Depth = vcat(CS_max_depth, CR_max_depth, SY_max_depth),
    Network = vcat(repeat(["CS"], length(CS_elapsed_time)), repeat(["CR"], length(CR_elapsed_time)), repeat(["SY"], length(SY_elapsed_time))),
)

CSV.write("regression.csv", df_reg)