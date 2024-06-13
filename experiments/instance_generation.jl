using PulseAlgorithm
const PA = PulseAlgorithm
using Distributions
include("util.jl")
include("modified_data_loader.jl")

# Set the seed for reproducibility
Random.seed!(1234)

# Fix CV and ρ for all experiments
ρ = 1.0 #this is fixed for our experiments
CV = 0.8 #this is fixed for our experiments
α = 0.9
γ = 0.4
max_depth = 2

### Chicago Regional ###

# Instance generation
network_name = "CR"
toll_factor_ChicagoRegional = 0.1 #taken from TransportationsNetworks
distance_factor_ChicagoRegional = 0.25 #taken from TransportationsNetworks

folder_path = raw"data\ChicagoRegional\shortest_paths"
graph = PA.load_graph_from_ta(raw"data\ChicagoRegional\ChicagoRegional_net.tntp", raw"data\ChicagoRegional\ChicagoRegional_flow.tntp",  "CR", CV, toll_factor_ChicagoRegional, distance_factor_ChicagoRegional)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

n = 10 # sample n nodes from graph.nodes uniformly
keys_list = collect(keys(graph.nodes))  
sampled_keys = sample(keys_list, n, replace=false)  

for key in sampled_keys
    write_shortest_paths(graph, graph.nodes[key].name, folder_path, "CR")
end


### Chicago Sketch ###
network_name = "CS"
toll_factor_ChicagoSketch = 0.02 #taken from TransportationsNetworks
distance_factor_ChicagoSketch = 0.04 #taken from TransportationsNetworks

folder_path = raw"data\ChicagoSketch\shortest_paths"
graph = PulseAlgorithm.load_graph_from_ta(raw"data\ChicagoSketch\ChicagoSketch_net.tntp", raw"data\ChicagoSketch\ChicagoSketch_flow.tntp", "CS", CV, toll_factor_ChicagoSketch, distance_factor_ChicagoSketch)
covariance_dict = get_covariance_dict(graph, ρ, max_depth)

n = 10 # sample n nodes from graph.nodes uniformly
keys_list = collect(keys(graph.nodes))
sampled_keys = sample(keys_list, n, replace=false)
println(sampled_keys)

for key in sampled_keys
    write_shortest_paths(graph, graph.nodes[key].name, folder_path, "CS")
end

### Austin ###
network_name = "AU"
toll_factor_Austin = 0.1 #Use Chicago Regional's toll factor as a reference
distance_factor_Austin = 0.25 #Use Chicago Regional's distance factor as a reference

folder_path = raw"data\Austin\shortest_paths"
graph = load_graph_from_ta_without_flow(raw"data\Austin\Austin_net.tntp", "AU", CV, toll_factor_Austin, distance_factor_Austin)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

n = 10 # sample n nodes from graph.nodes uniformly
keys_list = collect(keys(graph.nodes))
sampled_keys = sample(keys_list, n, replace=false)
println(sampled_keys)

for key in sampled_keys
    write_shortest_paths(graph, graph.nodes[key].name, folder_path, "AU")
end

### Sydney ###
network_name = "SY"
toll_factor_Sydney = 0.1 #Use Chicago Regional's toll factor as a reference
distance_factor_Sydney = 0.25 #Use Chicago Regional's distance factor as a reference

folder_path = raw"data\Sydney\shortest_paths"
graph = load_graph_from_ta_without_flow_Sydney(raw"data\Sydney\Sydney_net.tntp", "SY", CV, toll_factor_Sydney, distance_factor_Sydney)
covariance_dict = PA.get_covariance_dict(graph, ρ, max_depth)

n = 10 # sample n nodes from graph.nodes uniformly
keys_list = collect(keys(graph.nodes))
sampled_keys = sample(keys_list, n, replace=false)
println(sampled_keys)

for i in ProgressBar(1:length(sampled_keys))
    key = sampled_keys[i]
    write_shortest_paths(graph, graph.nodes[key].name, folder_path, "SY")
end