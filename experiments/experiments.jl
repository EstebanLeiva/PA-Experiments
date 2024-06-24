using PulseAlgorithm
const PA = PulseAlgorithm
using Distributions
using PlotlyJS, DataFrames
include("util.jl")
include("modified_data_loader.jl")

# Fixed parameters for our experiments
ρ = 1.0 
CV = 0.8 
α = 0.9
γ = 0.4
max_depth = 2

### Chicago Regional ###
network_name = "CR"
toll_factor_ChicagoRegional = 0.1 #taken from TransportationsNetworks
distance_factor_ChicagoRegional = 0.25 #taken from TransportationsNetworks

folder_path = raw"data\ChicagoRegional\shortest_paths"
graph = PA.load_graph_from_ta(raw"data\ChicagoRegional\ChicagoRegional_net.tntp", 
raw"data\ChicagoRegional\ChicagoRegional_flow.tntp",  "CR", CV, toll_factor_ChicagoRegional, distance_factor_ChicagoRegional)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

sampled_keys_CR = [10584, 9053, 5332, 10172, 6784, 2611, 8851, 11851, 8392, 11012]

n = 100

total_instance_info = aggregate_experiment_results(sampled_keys_CR, graph, ρ, α, γ, max_depth, folder_path, network_name, n)

avg_nondominated_paths = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys)*n)
avg_elapsed_time = total_instance_info["total_elapsed_time"] / (length(sampled_keys)*n)
avg_pruned_by_bounds = total_instance_info["pruned_by_bounds"] / (length(sampled_keys)*n)
avg_pruned_by_feasibility = total_instance_info["pruned_by_feasibility"] / (length(sampled_keys)*n)
avg_length_pruned_by_bounds = (total_instance_info["total_length_pruned_by_bounds"]/total_instance_info["pruned_by_bounds"])
avg_length_pruned_by_feasibility = (total_instance_info["total_length_pruned_by_feasibility"]/total_instance_info["pruned_by_feasibility"])

### Austin ###
network_name = "AU"
toll_factor_Austin = 0.1 #Use Chicago Regional's toll factor as a reference
distance_factor_Austin = 0.25 #Use Chicago Regional's distance factor as a reference

folder_path = raw"data\Austin\shortest_paths"
graph = load_graph_from_ta_without_flow(raw"data\Austin\Austin_net.tntp", "AU", CV, toll_factor_Austin, distance_factor_Austin)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

sampled_keys_AU = [149, 165, 252, 1220, 1573, 2561, 4529, 5389, 6192, 7352]

n = 100

total_instance_info = aggregate_experiment_results(sampled_keys_AU, graph, ρ, α, γ, max_depth, folder_path, network_name, n)

avg_nondominated_paths = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys)*n)
avg_elapsed_time = total_instance_info["total_elapsed_time"] / (length(sampled_keys)*n)
avg_pruned_by_bounds = total_instance_info["pruned_by_bounds"] / (length(sampled_keys)*n)
avg_pruned_by_feasibility = total_instance_info["pruned_by_feasibility"] / (length(sampled_keys)*n)
avg_length_pruned_by_bounds = (total_instance_info["total_length_pruned_by_bounds"]/total_instance_info["pruned_by_bounds"])
avg_length_pruned_by_feasibility = (total_instance_info["total_length_pruned_by_feasibility"]/total_instance_info["pruned_by_feasibility"])

### Sydney ###
network_name = "SY"
toll_factor_Sydney = 0.1 #Use Chicago Regional's toll factor as a reference
distance_factor_Sydney = 0.25 #Use Chicago Regional's distance factor as a reference

folder_path = raw"data\Sydney\shortest_paths"
graph = load_graph_from_ta_without_flow_Sydney(raw"data\Sydney\Sydney_net.tntp", "SY", CV, toll_factor_Sydney, distance_factor_Sydney)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

sampled_keys_SY = [3801, 8042, 11477, 13645, 15010, 15665, 16588, 16774, 17189, 17940]

n = 10

total_instance_info = aggregate_experiment_results(sampled_keys_SY, graph, ρ, α, γ, max_depth, folder_path, network_name, n)

avg_nondominated_paths = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys)*n)
avg_elapsed_time = total_instance_info["total_elapsed_time"] / (length(sampled_keys)*n)
avg_pruned_by_bounds = total_instance_info["pruned_by_bounds"] / (length(sampled_keys)*n)
avg_pruned_by_feasibility = total_instance_info["pruned_by_feasibility"] / (length(sampled_keys)*n)
avg_length_pruned_by_bounds = (total_instance_info["total_length_pruned_by_bounds"]/total_instance_info["pruned_by_bounds"])
avg_length_pruned_by_feasibility = (total_instance_info["total_length_pruned_by_feasibility"]/total_instance_info["pruned_by_feasibility"])