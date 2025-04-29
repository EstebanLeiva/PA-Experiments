mutable struct TA_Data
    network_name::String

    number_of_zones::Int
    number_of_nodes::Int
    first_thru_node::Int
    number_of_links::Int

    start_node::Array{Int,1}
    end_node::Array{Int,1}
    capacity::Array{Float64,1}
    link_length::Array{Float64,1}
    free_flow_time::Array{Float64,1}
    B::Array{Float64,1}
    power::Array{Float64,1}
    speed_limit::Array{Float64,1}
    toll::Array{Float64,1}
    link_type::Array{Int64,1}
end

function load_ta_network(network_data_file, network_name)
    search_sc(s, c) = something(findfirst(isequal(c), s), 0)

    @assert ispath(network_data_file)

    ##################################################
    # Network Data
    ##################################################


    number_of_zones = 0
    number_of_links = 0
    number_of_nodes = 0
    first_thru_node = 0

    n = open(network_data_file, "r")

    while (line = readline(n)) != ""
        if occursin("<NUMBER OF ZONES>", line)
            number_of_zones = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<NUMBER OF NODES>", line)
            number_of_nodes = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<FIRST THRU NODE>", line)
            first_thru_node = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<NUMBER OF LINKS>", line)
            number_of_links = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<END OF METADATA>", line)
            break
        end
    end

    @assert number_of_links > 0

    init_node = Array{Int64}(undef, number_of_links)
    term_node = Array{Int64}(undef, number_of_links)
    capacity = zeros(number_of_links)
    link_length = zeros(number_of_links)
    free_flow_time = zeros(number_of_links)
    b = zeros(number_of_links)
    power = zeros(number_of_links)
    speed_limit = zeros(number_of_links)
    toll = zeros(number_of_links)
    link_type = Array{Int64}(undef, number_of_links)

    idx = 1
    while !eof(n)
        line = readline(n)
        if occursin("~", line) || line == ""
            continue
        end

        if occursin(";", line)
            line = strip(line, [' ', '\n', ';'])
            line = replace(line, ";" => "")

            numbers = split(line)
            init_node[idx] = parse(Int64, numbers[1])
            term_node[idx] = parse(Int64, numbers[2])
            capacity[idx] = parse(Float64, numbers[3])
            link_length[idx] = parse(Float64, numbers[4])
            free_flow_time[idx] = parse(Float64, numbers[5])
            b[idx] = parse(Float64, numbers[6])
            power[idx] = parse(Float64, numbers[7])
            speed_limit[idx] = parse(Float64, numbers[8])
            toll[idx] = parse(Float64, numbers[9])
            link_type[idx] = Int(round(parse(Float64, numbers[10])))

            idx = idx + 1
        end
    end

    # Preparing data to return
    ta_data = TA_Data(
        network_name,
        number_of_zones,
        number_of_nodes,
        first_thru_node,
        number_of_links,
        init_node,
        term_node,
        capacity,
        link_length,
        free_flow_time,
        b,
        power,
        speed_limit,
        toll,
        link_type)

    return ta_data

end

### OUR FUNCTIONS ###
function load_ta_network_Sydney(network_data_file, network_name)
    search_sc(s, c) = something(findfirst(isequal(c), s), 0)

    @assert ispath(network_data_file)

    ##################################################
    # Network Data
    ##################################################


    number_of_zones = 0
    number_of_links = 0
    number_of_nodes = 0
    first_thru_node = 0

    n = open(network_data_file, "r")

    while (line = readline(n)) != ""
        if occursin("<NUMBER OF ZONES>", line)
            number_of_zones = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<NUMBER OF NODES>", line)
            number_of_nodes = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<FIRST THRU NODE>", line)
            first_thru_node = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<NUMBER OF LINKS>", line)
            number_of_links = parse(Int, line[search_sc(line, '>')+1:end])
        elseif occursin("<END OF METADATA>", line)
            break
        end
    end

    @assert number_of_links > 0

    init_node = Array{Int64}(undef, number_of_links)
    term_node = Array{Int64}(undef, number_of_links)
    capacity = zeros(number_of_links)
    link_length = zeros(number_of_links)
    free_flow_time = zeros(number_of_links)
    b = zeros(number_of_links)
    power = zeros(number_of_links)
    speed_limit = zeros(number_of_links)
    toll = zeros(number_of_links)
    link_type = Array{Int64}(undef, number_of_links)

    idx = 1
    while !eof(n)
        line = readline(n)
        if occursin("~", line) || line == ""
            continue
        end

        if !occursin(";", line)
            line = strip(line, [' ', '\n', ';'])
            line = replace(line, ";" => "")
            numbers = split(line)
            init_node[idx] = parse(Int64, numbers[1])
            term_node[idx] = parse(Int64, numbers[2])
            capacity[idx] = parse(Float64, numbers[3])
            link_length[idx] = parse(Float64, numbers[4])
            free_flow_time[idx] = parse(Float64, numbers[5])
            b[idx] = parse(Float64, numbers[6])
            power[idx] = parse(Float64, numbers[7])

            idx = idx + 1
        end
    end

    # Preparing data to return
    ta_data = TA_Data(
        network_name,
        number_of_zones,
        number_of_nodes,
        first_thru_node,
        number_of_links,
        init_node,
        term_node,
        capacity,
        link_length,
        free_flow_time,
        b,
        power,
        speed_limit,
        toll,
        link_type)

    return ta_data

end

function calculate_avg_fft_coefficient(ta_data::TA_Data)
    sum = 0
    for i in 1:length(ta_data.free_flow_time)
        sum += ta_data.free_flow_time[i] / ta_data.link_length[i]
    end
    return sum / length(ta_data.free_flow_time)
end

function load_graph_from_ta_without_flow(tntp_file_dir::String, network_name::String, CV::Float64, toll_factor::Float64, length_factor::Float64)
    ta_data = load_ta_network(tntp_file_dir, network_name)
    new_graph = PA.Graph(Dict{Int, PA.Node}(), Dict{String, Int}(), Dict{String, DefaultDict{Tuple{Int, Int, Int, Int}, Float64}}())

    for i in 1:length(ta_data.start_node)
        PA.find_or_add_node!(new_graph, string(ta_data.start_node[i]))
    end

    cost_flow = 100 # Like there is no flow data, we assume to be 100 for every link
    avg_fft_coefficient = calculate_avg_fft_coefficient(ta_data)

    for i in 1:length(ta_data.start_node)
        start = string(ta_data.start_node[i])
        dst = string(ta_data.end_node[i])
        if ta_data.free_flow_time[i] == 0
            fft = ta_data.link_length[i] * avg_fft_coefficient #Replace zero entries with the average free flow time coefficient * link length
        else
            fft = ta_data.free_flow_time[i]
        end
        mean = fft * (1 + ta_data.B[i] * (cost_flow / ta_data.capacity[i])^ta_data.power[i])
        variance = CV * abs(mean - fft)
        cost = mean + toll_factor * ta_data.toll[i] + length_factor * ta_data.link_length[i]
        
        deterministic_info = Dict{String, Float64}("cost" => cost)
        random_info = Dict{String, Dict{String, Float64}}("time" => Dict("mean" => mean, "variance" => variance))
        PA.add_link!(new_graph, start, dst, deterministic_info, random_info)
    end
    return new_graph
end

function load_graph_from_ta_without_flow_Sydney(tntp_file_dir::String, network_name::String, CV::Float64, toll_factor::Float64, length_factor::Float64)
    ta_data = load_ta_network_Sydney(tntp_file_dir, network_name)
    new_graph = PA.Graph(Dict{Int, PA.Node}(), Dict{String, Int}(), Dict{String, DefaultDict{Tuple{Int, Int, Int, Int}, Float64}}())

    for i in 1:length(ta_data.start_node)
        PA.find_or_add_node!(new_graph, string(ta_data.start_node[i]))
    end

    cost_flow = 100 # Like there is no flow data, we assume to be 100 for every link
    avg_fft_coefficient = calculate_avg_fft_coefficient(ta_data)

    for i in 1:length(ta_data.start_node)
        start = string(ta_data.start_node[i])
        dst = string(ta_data.end_node[i])
        if ta_data.free_flow_time[i] == 0
            fft = ta_data.link_length[i] * avg_fft_coefficient #Replace zero entries with the average free flow time coefficient * link length
        else
            fft = ta_data.free_flow_time[i]
        end
        mean = fft * (1 + ta_data.B[i] * (cost_flow / ta_data.capacity[i])^ta_data.power[i])
        variance = mean * CV
        cost = mean + toll_factor * ta_data.toll[i] + length_factor * ta_data.link_length[i]
        
        deterministic_info = Dict{String, Float64}("cost" => cost)
        random_info = Dict{String, Dict{String, Float64}}("time" => Dict("mean" => mean, "variance" => variance))
        PA.add_link!(new_graph, start, dst, deterministic_info, random_info)
    end
    return new_graph
end