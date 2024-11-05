using PulseAlgorithm
using Distributions
using Random
using ProgressBars
using StatsPlots
using DataFrames
using CSV
const PA = PulseAlgorithm
include("util.jl")
include("sydney_loader.jl")

Random.seed!(1234)
path = raw"C:\Users\esteb\OneDrive - Universidad de los andes\Documentos\ANDES\Noveno Semestre\PA-RSPPs\PA-Experiments"

ρ = 1.0
CV = 0.8
α = 0.9
γ = 0.4
max_depth = 2
n = 100  #number of start_nodes per target_node

## Load networks ##

# Chicago Sketch (CS)
net_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_net.tntp")
flow_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_flow.tntp")

toll_factor_CS = 0.02
distance_factor_CS = 0.04

println(net_CS)
graph_CS = PA.load_graph_from_ta(net_CS, flow_CS, "CS", CV, toll_factor_CS, distance_factor_CS)
cov_CS = PA.get_covariance_dict(graph_CS, ρ, max_depth)

# Chicago Regional (CR)
net_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_net.tntp")
flow_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_flow.tntp")

toll_factor_CR = 0.1
distance_factor_CR = 0.25

graph_CR = PA.load_graph_from_ta(net_CR, flow_CR, "CR", CV, toll_factor_CR, distance_factor_CR)
cov_CR = PA.get_covariance_dict(graph_CR, ρ, max_depth)

# Sydney (SY)
net_SY = joinpath(path, raw"data\Sydney\sydney_net.tntp")
flow_SY = joinpath(path, raw"data\Sydney\sydney_flow.tntp")

toll_factor_SY = 0.1 #Use Chicago Regional's toll factor as a reference
distance_factor_SY = 0.25 #Use Chicago Regional's toll factor as a reference

graph_SY = load_graph_from_ta_without_flow_Sydney(net_SY, "SY", CV, toll_factor_SY, distance_factor_SY)
cov_SY = PA.get_covariance_dict(graph_SY, ρ, max_depth)

## Computational Time Experiments ##

# CS
sampled_keys_CS = sample(collect(keys(graph_CS.nodes)), 10, replace=false)
println("Sampled keys for Chicago Sketch: ", sampled_keys_CS)
info_CS = aggregate_experiments(sampled_keys_CS, graph_CS, α, γ, cov_CS, n)

CS_avg_nondominated_paths = info_CS["number_nondominanted_paths"]
CS_avg_elapsed_time = info_CS["total_elapsed_time"]
CS_avg_pruned_by_bounds = info_CS["pruned_by_bounds"]
CS_avg_pruned_by_feasibility = info_CS["pruned_by_feasibility"]
CS_avg_length_pruned_by_bounds = (info_CS["total_length_pruned_by_bounds"] ./ info_CS["pruned_by_bounds"])
CS_avg_length_pruned_by_feasibility = (info_CS["total_length_pruned_by_feasibility"] ./ info_CS["pruned_by_feasibility"])

println("---------------------------------")
println("Chicago Sketch")


# CR
sampled_keys_CR = sample(collect(keys(graph_CR.nodes)), 10, replace=false)
println("Sampled keys for Chicago Regional: ", sampled_keys_CR)
info_CR = aggregate_experiments(sampled_keys_CR, graph_CR, α, γ, cov_CR, n)

CR_avg_nondominated_paths = info_CR["number_nondominanted_paths"]
CR_avg_elapsed_time = info_CR["total_elapsed_time"]
CR_avg_pruned_by_bounds = info_CR["pruned_by_bounds"]
CR_avg_pruned_by_feasibility = info_CR["pruned_by_feasibility"]
CR_avg_length_pruned_by_bounds = (info_CR["total_length_pruned_by_bounds"] ./ info_CR["pruned_by_bounds"])
CR_avg_length_pruned_by_feasibility = (info_CR["total_length_pruned_by_feasibility"] ./ info_CR["pruned_by_feasibility"])

println("---------------------------------")
println("Chicago Regional")


# SY
sampled_keys_SY = sample(collect(keys(graph_SY.nodes)), 10, replace=false)
println("Sampled keys for Sydney: ", sampled_keys_SY)
info_SY = aggregate_experiments(sampled_keys_SY, graph_SY, α, γ, cov_SY, n)
SY_avg_nondominated_paths = info_SY["number_nondominanted_paths"]
SY_avg_elapsed_time = info_SY["total_elapsed_time"]
SY_avg_pruned_by_bounds = info_SY["pruned_by_bounds"]
SY_avg_pruned_by_feasibility = info_SY["pruned_by_feasibility"]
SY_avg_length_pruned_by_bounds = (info_SY["total_length_pruned_by_bounds"] ./ info_SY["pruned_by_bounds"])
SY_avg_length_pruned_by_feasibility = (info_SY["total_length_pruned_by_feasibility"] ./ info_SY["pruned_by_feasibility"])

println("---------------------------------")
println("Sydney")

println("---------------------------------")
println("Plotting")

### DFs ###

df_time = DataFrame(CS_t = CS_avg_elapsed_time, CR_t = CR_avg_elapsed_time, SY_b = SY_avg_elapsed_time)
CSV.write("time_sarp.csv", df_time)

df_prune = DataFrame(CS_b = CS_avg_pruned_by_bounds, 
                        CS_i = CS_avg_pruned_by_feasibility,
                        CR_b = CR_avg_pruned_by_bounds, 
                        CR_i = CR_avg_pruned_by_feasibility,
                        SY_b = SY_avg_pruned_by_bounds,
                        SY_i = SY_avg_pruned_by_feasibility)

CSV.write("prune_sarp.csv", df_prune)

df_length_prune = DataFrame(CS_b = CS_avg_length_pruned_by_bounds, 
                        CS_i = CS_avg_length_pruned_by_feasibility,
                        CR_b = CR_avg_length_pruned_by_bounds, 
                        CR_i = CR_avg_length_pruned_by_feasibility,
                        SY_b = SY_avg_length_pruned_by_bounds,
                        SY_i = SY_avg_length_pruned_by_feasibility)

CSV.write("length_prune_sarp.csv", df_length_prune)