using PulseAlgorithm
const PA = PulseAlgorithm
using Random
using Distributions
include("util_sdrspp.jl")
include("sydney_loader.jl")
include("chen_loader.jl")

Random.seed!(1234)
path = raw"C:\Users\investigacion\Documents\PA-Experiments"


α = 0.9
max_depth = 1
n = 10 #number of start_nodes per target_node

## Load networks ##

# Chicago Sketch (CS)
net_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_net.tntp")
flow_CS = joinpath(path, raw"data\ChicagoSketch\ChicagoSketch_flow.tntp")
folder_CS = joinpath(path, raw"data\ChicagoSketch")

divisor_CS = 5280.0
max_speed_CS = 100.0

graph_CS = load_graph_from_ta_chen(net_CS, "")
cov_CS = get_covariance_dict_chen(graph_CS, max_depth)

# Chicago Regional (CR)
net_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_net.tntp")
flow_CR = joinpath(path, raw"data\ChicagoRegional\ChicagoRegional_flow.tntp")
folder_CR = joinpath(path, raw"data\ChicagoRegional")

divisor_CR = 5280.0
max_speed_CR = 100.0

graph_CR = load_graph_from_ta_chen(net_CR, "")
cov_CR = get_covariance_dict_chen(graph_CR, max_depth)

# Sydney (SY)
net_SY = joinpath(path, raw"data\Sydney\sydney_net.tntp")
folder_SY = joinpath(path, raw"data\Sydney")

divisor_SY = 1.0
max_speed_SY = 100.0

graph_SY = load_graph_Sydney_chen(net_SY, "")
cov_SY = get_covariance_dict_chen(graph_SY, max_depth)

## Computational Time Experiments ##

# CS
sampled_keys_CS = sample(collect(keys(graph_CS.nodes)), 10, replace=false)
println("Sampled keys for Chicago Sketch: ", sampled_keys_CS)
pulse_info_CS, erspa_info_CS = aggregate_experiments(sampled_keys_CS, graph_CS, α, cov_CS, folder_CS, true, n, max_speed_CS, divisor_CS)

CS_avg_nondominated_paths = pulse_info_CS["number_nondominanted_paths"] / (length(sampled_keys_CS) * n)
CS_avg_elapsed_time = pulse_info_CS["total_elapsed_time"] / (length(sampled_keys_CS) * n)
CS_avg_pruned_by_bounds = pulse_info_CS["pruned_by_bounds"] / (length(sampled_keys_CS) * n)
CS_avg_length_pruned_by_bounds = (pulse_info_CS["total_length_pruned_by_bounds"] / pulse_info_CS["pruned_by_bounds"])
CS_avg_nondominated_paths_erspa = erspa_info_CS["number_nondominanted_paths"] / (length(sampled_keys_CS) * n)
CS_avg_elapsed_time_erspa = erspa_info_CS["total_elapsed_time"] / (length(sampled_keys_CS) * n)

println("-------------------------------------------------")
println("Chicago Sketch")
println("PULSE CS_avg_nondominated_paths: ", CS_avg_nondominated_paths)
println("PULSE CS_avg_elapsed_time: ", CS_avg_elapsed_time)
println("PULSE CS_avg_pruned_by_bounds: ", CS_avg_pruned_by_bounds)
println("PULSE CS_avg_length_pruned_by_bounds: ", CS_avg_length_pruned_by_bounds)
println("ERSPA CS_avg_nondominated_paths: ", CS_avg_nondominated_paths_erspa)
println("ERSPA CS_avg_elapsed_time: ", CS_avg_elapsed_time_erspa)

# CR

sampled_keys_CR = sample(collect(keys(graph_CR.nodes)), 10, replace=false)
println("Sampled keys for Chicago Regional: ", sampled_keys_CR)
pulse_info_CR, erspa_info_CR = aggregate_experiments(sampled_keys_CR, graph_CR, α, cov_CR, folder_CR, true, n, max_speed_CR, divisor_CR)

CR_avg_nondominated_paths = pulse_info_CR["number_nondominanted_paths"] / (length(sampled_keys_CR) * n)
CR_avg_elapsed_time = pulse_info_CR["total_elapsed_time"] / (length(sampled_keys_CR) * n)
CR_avg_pruned_by_bounds = pulse_info_CR["pruned_by_bounds"] / (length(sampled_keys_CR) * n)
CR_avg_length_pruned_by_bounds = (pulse_info_CR["total_length_pruned_by_bounds"] / pulse_info_CR["pruned_by_bounds"])
CR_avg_nondominated_paths_erspa = erspa_info_CR["number_nondominanted_paths"] / (length(sampled_keys_CR) * n)
CR_avg_elapsed_time_erspa = erspa_info_CR["total_elapsed_time"] / (length(sampled_keys_CR) * n)

println("-------------------------------------------------")
println("Chicago Regional")
println("PULSE CR_avg_nondominated_paths: ", CR_avg_nondominated_paths)
println("PULSE CR_avg_elapsed_time: ", CR_avg_elapsed_time) 
println("PULSE CR_avg_pruned_by_bounds: ", CR_avg_pruned_by_bounds)
println("PULSE CR_avg_length_pruned_by_bounds: ", CR_avg_length_pruned_by_bounds)
println("ERSPA CR_avg_nondominated_paths: ", CR_avg_nondominated_paths_erspa)
println("ERSPA CR_avg_elapsed_time: ", CR_avg_elapsed_time_erspa)

# SY
sampled_keys_SY = sample(collect(keys(graph_SY.nodes)), 10, replace=false)
println("Sampled keys for Sydney: ", sampled_keys_SY)
pulse_info_SY, erspa_info_SY = aggregate_experiments(sampled_keys_SY, graph_SY, α, cov_SY, folder_SY, true, n, max_speed_SY, divisor_SY)

SY_avg_nondominated_paths = pulse_info_SY["number_nondominanted_paths"] / (length(sampled_keys_SY) * n)
SY_avg_elapsed_time = pulse_info_SY["total_elapsed_time"] / (length(sampled_keys_SY) * n)
SY_avg_pruned_by_bounds = pulse_info_SY["pruned_by_bounds"] / (length(sampled_keys_SY) * n)
SY_avg_length_pruned_by_bounds = (pulse_info_SY["total_length_pruned_by_bounds"] / pulse_info_SY["pruned_by_bounds"])
SY_avg_nondominated_paths_erspa = erspa_info_SY["number_nondominanted_paths"] / (length(sampled_keys_SY) * n)
SY_avg_elapsed_time_erspa = erspa_info_SY["total_elapsed_time"] / (length(sampled_keys_SY) * n)

println("-------------------------------------------------")
println("Sydney")
println("PULSE SY_avg_nondominated_paths: ", SY_avg_nondominated_paths)
println("PULSE SY_avg_elapsed_time: ", SY_avg_elapsed_time)
println("PULSE SY_avg_pruned_by_bounds: ", SY_avg_pruned_by_bounds)
println("PULSE SY_avg_length_pruned_by_bounds: ", SY_avg_length_pruned_by_bounds)
println("ERSPA SY_avg_nondominated_paths: ", SY_avg_nondominated_paths_erspa)
println("ERSPA SY_avg_elapsed_time: ", SY_avg_elapsed_time_erspa)