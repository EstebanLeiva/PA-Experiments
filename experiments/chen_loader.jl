include("sydney_loader.jl")
include("util.jl")
function load_graph_from_ta_chen(tntp_file_dir::String, network_name::String)
    ta_data = PA.load_ta(tntp_file_dir, network_name)
    new_graph = PA.Graph(Dict{Int, PA.Node}(), Dict{String, Int}(), Dict{String, DefaultDict{Tuple{Int, Int, Int, Int}, Float64}}())

    for i in 1:length(ta_data.start_node)
        PA.find_or_add_node!(new_graph, string(ta_data.start_node[i]))
    end

    for i in 1:length(ta_data.start_node)
        start = string(ta_data.start_node[i])
        dst = string(ta_data.end_node[i])
        length = ta_data.link_length[i]
        speed = rand(Uniform(10, 100))
        mean = length / speed
        CV = rand(Uniform(0.1, 1))
        variance = mean * CV
        cost = 0.0
        deterministic_info = Dict{String, Float64}("cost" => cost)
        random_info = Dict{String, Dict{String, Float64}}("time" => Dict("mean" => mean, "variance" => variance))
        PA.add_link!(new_graph, start, dst, deterministic_info, random_info)
    end
    return new_graph
end

function get_covariance_dict_chen(graph::PA.Graph, random_variable::String, max_depth::Int)
    covariance_dict = Dict{String, DefaultDict{Tuple{Int, Int, Int, Int}, Float64}}() 
    covariance_dict[random_variable] = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0)
    links = PA.get_links_info(graph)
    for link in keys(links)
        visited_pairlinks = modified_dfs(graph, link, max_depth, 1, Dict{Tuple{Int, Int}, Int}(), -1)
        for pairlink in keys(visited_pairlinks)
            if pairlink != link 
                covariance_dict[random_variable][(link[1], link[2], pairlink[1], pairlink[2])] = rand(Uniform(0.1, 1)) * √(links[(pairlink)][2][random_variable]["variance"]) * √(links[link][2][random_variable]["variance"]) 
            end
            if pairlink == link
                covariance_dict[random_variable][(link[1], link[2], pairlink[1], pairlink[2])] = 1.0
            end
        end
    end
    return covariance_dict
end

function load_graph_Sydney_chen(tntp_file_dir::String, network_name::String)
    ta_data = load_ta_network_Sydney(tntp_file_dir, network_name)
    new_graph = PA.Graph(Dict{Int, PA.Node}(), Dict{String, Int}(), Dict{String, DefaultDict{Tuple{Int, Int, Int, Int}, Float64}}())

    for i in 1:length(ta_data.start_node)
        PA.find_or_add_node!(new_graph, string(ta_data.start_node[i]))
    end

    for i in 1:length(ta_data.start_node)
        start = string(ta_data.start_node[i])
        dst = string(ta_data.end_node[i])
        mean = rand(Uniform(10, 100))
        CV = rand(Uniform(0.1, 1))
        variance = mean * CV
        cost = 0.0
        deterministic_info = Dict{String, Float64}("cost" => cost)
        random_info = Dict{String, Dict{String, Float64}}("time" => Dict("mean" => mean, "variance" => variance))
        PA.add_link!(new_graph, start, dst, deterministic_info, random_info)
    end
    return new_graph
end

function load_graph_Philadelphia_chen(tntp_file_dir::String, network_name::String)
    ta_data = PA.load_ta(tntp_file_dir, network_name)
    new_graph = PA.Graph(Dict{Int, PA.Node}(), Dict{String, Int}(), Dict{String, DefaultDict{Tuple{Int, Int, Int, Int}, Float64}}())

    for i in 1:length(ta_data.start_node)
        PA.find_or_add_node!(new_graph, string(ta_data.start_node[i]))
    end 
    
    for i in 1:length(ta_data.start_node)
        start = string(ta_data.start_node[i])
        dst = string(ta_data.end_node[i])
        mean = rand(Uniform(10, 100))
        CV = rand(Uniform(0.1, 1))
        variance = mean * CV
        cost = 0.0
        deterministic_info = Dict{String, Float64}("cost" => cost)
        random_info = Dict{String, Dict{String, Float64}}("time" => Dict("mean" => mean, "variance" => variance))
        PA.add_link!(new_graph, start, dst, deterministic_info, random_info)
    end
    return new_graph
end