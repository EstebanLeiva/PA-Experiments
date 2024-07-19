using PulseAlgorithm
const PA = PulseAlgorithm
using Distributions
using PlotlyJS, DataFrames
include("util.jl")
include("util_negcorr.jl")
include("modified_data_loader.jl")

# Fixed parameters for our experiments
ρ = 1.0 
CV = 0.8 
α = 0.9
γ = 0.4
max_depth = 4

### Chicago Sketch ###
network_name = "CS"
toll_factor_ChicagoSketch = 0.02 #taken from TransportationsNetworks
distance_factor_ChicagoSketch = 0.04 #taken from TransportationsNetworks

folder_path = raw"data\ChicagoSketch\shortest_paths"
graph = PA.load_graph_from_ta(raw"data\ChicagoSketch\ChicagoSketch_net.tntp", raw"data\ChicagoSketch\ChicagoSketch_flow.tntp", "CS", CV, toll_factor_ChicagoSketch, distance_factor_ChicagoSketch)
covariance_dict = get_neg_covariance_dict(graph, ρ, max_depth)

sampled_keys_CS = [39, 195, 204, 440, 493, 682, 706, 767, 798, 812]

n = 1000

total_instance_info = aggregate_neg_experiment_results(sampled_keys_CS, graph, ρ, α, γ, max_depth, folder_path, network_name, n)

CS_avg_reliability = total_instance_info["post_reliability"] / (length(sampled_keys_CS)*n)


### Chicago Regional ###
network_name = "CR"
toll_factor_ChicagoRegional = 0.1 #taken from TransportationsNetworks
distance_factor_ChicagoRegional = 0.25 #taken from TransportationsNetworks

folder_path = raw"data\ChicagoRegional\shortest_paths"
graph = PA.load_graph_from_ta(raw"data\ChicagoRegional\ChicagoRegional_net.tntp",
raw"data\ChicagoRegional\ChicagoRegional_flow.tntp",  "CR", CV, toll_factor_ChicagoRegional, distance_factor_ChicagoRegional)
covariance_dict = get_neg_covariance_dict(graph, ρ, max_depth)

sampled_keys_CR = [10584, 9053, 5332, 10172, 6784, 2611, 8851, 11851, 8392, 11012]

n = 1000

total_instance_info = aggregate_neg_experiment_results(sampled_keys_CR, graph, ρ, α, γ, max_depth, folder_path, network_name, n)

CR_avg_reliability = total_instance_info["post_reliability"] / (length(sampled_keys_CR)*n)
