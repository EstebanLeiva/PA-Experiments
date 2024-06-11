using PulseAlgorithm
using Distributions

include("util.jl")

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
graph = PulseAlgorithm.load_graph_from_ta(raw"data\ChicagoRegional\ChicagoRegional_net.tntp", 
raw"data\ChicagoRegional\ChicagoRegional_flow.tntp",  "CR", CV, toll_factor_ChicagoRegional, distance_factor_ChicagoRegional)
covariance_dict = get_covariance_dict(graph, ρ, max_depth)

n = 10 # sample n nodes from graph.nodes uniformly
keys_list = collect(keys(graph.nodes))  
sampled_keys = sample(keys_list, n, replace=false)  

for key in sampled_keys
    PulseAlgorithm.write_shortest_paths(graph, graph.nodes[key].name, folder_path, "CR")
end