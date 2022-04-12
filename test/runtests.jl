using Test
using JSON
using CSV
using DataFrames
include("../src/GHPGHX.jl")

inputs_dict = JSON.parsefile("inputs/inputs_base.json")
cop_map_df = CSV.read("inputs/cop_map.csv", DataFrame)
# Generate a "records" style dictionary from the 
cop_map_list = []
for i in eachrow(cop_map_df)
    dict_record = Dict(names(i)[k]=>i[names(i)[k]] for k in 1:length(names(i)))
    push!(cop_map_list, dict_record)
end
inputs_dict["cop_map_eft_heating_cooling"] = cop_map_list

# Heating and cooling loads
inputs_dict["heating_thermal_load_mmbtu_per_hr"] = ones(8760) * 2
inputs_dict["cooling_thermal_load_ton"] = ones(8760) * 300

# TODO Add PVWatts API call for dry bulb temperature based on lat/long
inputs_dict["ambient_temperature_f"] = ones(8760) * 59.0

# TODO Add ground_k_by_climate_zone lookup functionality from REopt.jl and need lat/long too
inputs_dict["ground_thermal_conductivity_btu_per_hr_ft_f"] = 1.2

@info "Starting GHPGHX" #with timeout of $(timeout) seconds..."
results, inputs_params = GHPGHX.ghp_model(inputs_dict)
# Create a dictionary of the results data needed for REopt
ghpghx_results = GHPGHX.get_ghpghx_results_for_reopt(results, inputs_params)
@info "GHPGHX model solved" #with status $(results["status"])."


#@test 1 == 1