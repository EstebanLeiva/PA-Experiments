using PulseAlgorithm
const PA = PulseAlgorithm
using Distributions
using PlotlyJS, DataFrames
include("util_sdrspp.jl")
include("modified_data_loader.jl")

# Fixed parameters for our experiments
ρ = 1.0 
CV = 0.8 
α = 0.9
max_depth = 1

network_name = "CS"
toll_factor_ChicagoSketch = 0.02 #taken from TransportationsNetworks
distance_factor_ChicagoSketch = 0.04 #taken from TransportationsNetworks

folder_path = raw"data\ChicagoSketch\shortest_paths"
graph = PA.load_graph_from_ta(raw"data\ChicagoSketch\ChicagoSketch_net.tntp", raw"data\ChicagoSketch\ChicagoSketch_flow.tntp", "CS", CV, toll_factor_ChicagoSketch, distance_factor_ChicagoSketch)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

sampled_keys_CS = [39, 195, 204, 440, 493, 682, 706, 767, 798, 812]

n = 1000

total_instance_info = aggregate_experiment_results(sampled_keys_CS, graph, ρ, α, max_depth, folder_path, network_name, n)

CS_avg_nondominated_paths = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys_CS)*n)
CS_avg_elapsed_time = total_instance_info["total_elapsed_time"] / (length(sampled_keys_CS)*n)
CS_avg_pruned_by_bounds = total_instance_info["pruned_by_bounds"] / (length(sampled_keys_CS)*n)
CS_avg_length_pruned_by_bounds = (total_instance_info["total_length_pruned_by_bounds"]/total_instance_info["pruned_by_bounds"])

total_instance_info = erspa_aggregate_experiment_results(sampled_keys_CS, graph, ρ, α, max_depth, folder_path, network_name, n)

CS_avg_nondominated_paths_erspa = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys_CS)*n)
CS_avg_elapsed_time_erspa = total_instance_info["total_elapsed_time"] / (length(sampled_keys_CS)*n)

### Chicago Regional ###
network_name = "CR"
toll_factor_ChicagoRegional = 0.1 #taken from TransportationsNetworks
distance_factor_ChicagoRegional = 0.25 #taken from TransportationsNetworks

folder_path = raw"data\ChicagoRegional\shortest_paths"
graph = PA.load_graph_from_ta(raw"data\ChicagoRegional\ChicagoRegional_net.tntp", raw"data\ChicagoRegional\ChicagoRegional_flow.tntp",  "CR", CV, toll_factor_ChicagoRegional, distance_factor_ChicagoRegional)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

sampled_keys_CR = [10584, 9053, 5332, 10172, 6784, 2611, 8851, 11851, 8392, 11012]

n = 1

total_instance_info = aggregate_experiment_results(sampled_keys_CR, graph, ρ, α, max_depth, folder_path, network_name, n)

CR_avg_nondominated_paths = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys_CR)*n)
CR_avg_elapsed_time = total_instance_info["total_elapsed_time"] / (length(sampled_keys_CR)*n)
CR_avg_pruned_by_bounds = total_instance_info["pruned_by_bounds"] / (length(sampled_keys_CR)*n)
CR_avg_length_pruned_by_bounds = (total_instance_info["total_length_pruned_by_bounds"]/total_instance_info["pruned_by_bounds"])

total_instance_info = erspa_aggregate_experiment_results(sampled_keys_CR, graph, ρ, α, max_depth, folder_path, network_name, n)

CR_avg_nondominated_paths_erspa = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys_CR)*n)
CR_avg_elapsed_time_erspa = total_instance_info["total_elapsed_time"] / (length(sampled_keys_CR)*n)
