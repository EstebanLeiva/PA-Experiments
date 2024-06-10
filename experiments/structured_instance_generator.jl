using PulseAlgorithm
using DataStructures

function modified_dfs(graph::Graph, start_link::Tuple{Int, Int}, max_depth::Int, depth::Int, visited_pairlinks::Dict{Tuple{Int, Int}, Int}, previous_node::Int)
    if depth > max_depth
        return visited_pairlinks
    end

    for adjacent in keys(graph.nodes[start_link[2]].links)
        if previous_node != adjacent
            if haskey(visited_pairlinks, (start_link[2], adjacent))
                visited_pairlinks[(start_link[2], adjacent)] = min(visited_pairlinks[(start_link[2], adjacent)], depth)
            else
                visited_pairlinks[(start_link[2], adjacent)] = depth
            end
            modified_dfs(graph, (start_link[2], adjacent), max_depth, depth + 1, visited_pairlinks, adjacent)
        end
    end       

    return visited_pairlinks
end

function get_covariance_dict(graph::Graph, ρ::Float64, max_depth::Int)
    covariance_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0) #default value of dic is 0.0
    links = get_links_info(graph)
    for link in keys(links)
        visited_pairlinks = modified_dfs(graph, link, max_depth, 1, Dict{Tuple{Int, Int}, Int}(), -1)
        for pairlink in keys(visited_pairlinks)
            if pairlink != link 
                covariance_dict[(link[1], link[2], pairlink[1], pairlink[2])] = (ρ/(visited_pairlinks[pairlink])) * √(links[(pairlink)][3]) * √(links[link][3]) #Corredor structure of ρ's
            end
            if pairlink == link
                covariance_dict[(link[1],link[2],pairlink[1],pairlink[2])] = 1.0
            end
        end
    end
    return covariance_dict
end

function get_quantile_path(graph::Graph, path::Vector{Int}, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
    mean = 0.0
    variance = 0.0
    covariance_term = 0.0
    for i in 1:length(path)-1
        mean += graph.nodes[path[i]].links[path[i+1]].mean
        variance += graph.nodes[path[i]].links[path[i+1]].variance
        for ii in i + 1:length(path)-1
            covariance_term += 2*cov_dict[(path[i], path[i+1], path[ii], path[ii+1])]
        end
    end
    return mean, variance, covariance_term
end

function get_timeBudget(graph::Graph, start_node::Int, target_node::Int, α::Float64, γ::Float64, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
    shortest_mean_path = dijkstra_between2Nodes(graph, start_node, target_node, "mean")
    shortest_cost_path = dijkstra_between2Nodes(graph, start_node, target_node, "cost")
    
    cost_min_mean = 0.0
    for i in 1:length(shortest_mean_path)-1
        cost_min_mean = cost_min_mean + graph.nodes[shortest_mean_path[i]].links[shortest_mean_path[i+1]].cost 
    end
    
    mean, variance, covariance_term = get_quantile_path(graph, shortest_mean_path, cov_dict)
    T_t_α = quantile(Normal(mean, √(variance+covariance_term)), α)

    mean, variance, covariance_term = get_quantile_path(graph, shortest_cost_path, cov_dict)
    T_c_α = quantile(Normal(mean, √(variance+covariance_term)), α)

    T = T_t_α + (T_c_α - T_t_α) * (1 - γ)
    return T, shortest_mean_path, cost_min_mean
end

function write_shortest_paths(graph::Graph, target_node::String, folder_path::String, network_name::String)
    variance_costs = dijkstra(graph, target_node, "variance")
    file_path = joinpath(folder_path, "variance_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(variance_costs = variance_costs), writeheader = false)

    mean_costs = dijkstra(graph, target_node, "mean")
    file_path = joinpath(folder_path, "mean_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(mean_costs = mean_costs), writeheader = false)

    minimum_costs = dijkstra(graph, target_node, "cost")
    file_path = joinpath(folder_path, "minimum_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(minimum_costs = minimum_costs), writeheader = false)
end

function preprocess_experiments(sp::SPulseGraph, folder_path::String, network_name::String, target_node::String)
    file_path = joinpath(folder_path, "minimum_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.minimum_costs =  data[:, 1] |> collect

    file_path = joinpath(folder_path, "mean_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.mean_costs =  data[:, 1] |> collect

    file_path = joinpath(folder_path, "variance_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.variance_costs =  data[:, 1] |> collect

    nodes = sort(collect(sp.G.nodes), by=x->sp.minimum_costs[x[1]], rev = true)
    possible_start_nodes = filter(x->sp.minimum_costs[x[1]] != Inf, nodes)

    sp.source_node = possible_start_nodes[1][2].name #the node farthest from target

    return sp.source_node
end

function run_experiments_time(graph::Graph, source_node::String, target_node::String, ρ::Float64, α::Float64, γ::Float64, max_depth::Int, folder_path::String, network_name::String, initial_bound::Bool)
    covariance_dict = get_covariance_dict(graph, ρ, max_depth)
    pulse = create_SPulseGraph(graph, α, covariance_dict, source_node, target_node, 0.0)
    source_node = preprocess_experiments(pulse, folder_path, network_name, target_node) 
    T, shortest_mean_path, cost_min_mean = get_timeBudget(graph, pulse.G.name_to_index[source_node], pulse.G.name_to_index[target_node], α, γ, covariance_dict)
    pulse.T_max = T
    
    if prob >= α && initial_bound
        elapsed_time = @elapsed begin
            run_pulse(pulse, shortest_mean_path, cost_min_mean)
        end
    else
        elapsed_time = @elapsed begin
            run_pulse(pulse)
        end
    end
    return elapsed_time, pulse.instance_info, (source_node, target_node), pulse.T_max, pulse.optimal_path
end