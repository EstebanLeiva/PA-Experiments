using PulseAlgorithm
const PA = PulseAlgorithm
using Distributions
using PlotlyJS, DataFrames
include("util.jl")

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
length(graph.nodes)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

sampled_keys = [10584, 9053, 5332, 10172, 6784, 2611, 8851, 11851, 8392, 11012]

n = 100

total_instance_info = aggregate_experiment_results(sampled_keys, graph, ρ, α, γ, max_depth, folder_path, network_name, n)

avg_nondominated_paths = total_instance_info["number_nondominanted_paths"] / (length(sampled_keys)*n)
avg_elapsed_time = total_instance_info["total_elapsed_time"] / (length(sampled_keys)*n)
avg_pruned_by_bounds = total_instance_info["pruned_by_bounds"] / (length(sampled_keys)*n)
avg_pruned_by_feasibility = total_instance_info["pruned_by_feasibility"] / (length(sampled_keys)*n)
avg_length_pruned_by_bounds = (total_instance_info["total_length_pruned_by_bounds"]/total_instance_info["pruned_by_bounds"])
avg_length_pruned_by_feasibility = (total_instance_info["total_length_pruned_by_feasibility"]/total_instance_info["pruned_by_feasibility"])

