include("util_table.jl")
Random.seed!(123)
path = raw"C:\Users\esteb\OneDrive - Universidad de los andes\Documentos\ANDES\Noveno Semestre\PA-RSPPs\PA-Experiments"

ρ = 1.0
CV = 0.8
α = 0.9
γ = 0.4
max_depth = 2

experiment_table(α, ρ, γ, max_depth)