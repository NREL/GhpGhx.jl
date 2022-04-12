using Test
using JSON
using CSV
include("../src/GHPGHX.jl")

inputs_dict = JSON.parsefile("test/inputs/inputs_base.json")
cop_map_df = CSV.File("test/inputs/cop_map.csv")
@info "Starting GHPGHX" #with timeout of $(timeout) seconds..."
results, inputs_params = ghp_model(inputs_dict)
# Create a dictionary of the results data needed for REopt
ghpghx_results = get_ghpghx_results_for_reopt(results, inputs_params)
@info "GHPGHX model solved" #with status $(results["status"])."


@test 1 == 1