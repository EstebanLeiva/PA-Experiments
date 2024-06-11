using PulseAlgorithm
using Distributions

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
graph = PulseAlgorithm.load_graph_from_ta(raw"data\ChicagoRegional\ChicagoRegional_net.tntp", 
raw"data\ChicagoRegional\ChicagoRegional_flow.tntp",  "CR", CV, toll_factor_ChicagoRegional, distance_factor_ChicagoRegional)
covariance_dict = PulseAlgorithm.get_covariance_dict(graph, ρ, max_depth)

sampled_keys = [10584, 9053, 5332, 10172, 6784, 2611, 8851, 11851, 8392, 11012]
target_node = graph.nodes[sampled_keys[1]].name

n = 10000
total_instance_info, total_elapsed_time = run_aggregated_experiments(graph, target_node, ρ, α, γ, max_depth, folder_path, network_name, true, n)

avg_elapsed_time = total_elapsed_time/n
#avg total instance instance_info
total_instance_info["pruned_by_bounds"]/n
total_instance_info["pruned_by_feasibility"]/n
total_instance_info["total_length_pruned_by_bounds"]/n
total_instance_info["total_length_pruned_by_feasibility"]/n