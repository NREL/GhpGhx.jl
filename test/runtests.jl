using Test
using JSON
using CSV
using DataFrames
include("../src/GhpGhx.jl")

@testset "WSHP system sizing" begin

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

    @info "Starting GhpGhx" #with timeout of $(timeout) seconds..."
    results, inputs_params = GhpGhx.ghp_model(inputs_dict)

    # Create a dictionary of the results data needed for REopt
    GhpGhx_results = GhpGhx.get_results_for_reopt(results, inputs_params)
    @info "GhpGhx model solved" #with status $(results["status"])."
 
    @test typeof(results) == GhpGhx.ResultsStruct
    @test GhpGhx_results["number_of_boreholes"] ≈ 1330 rtol=0.02
    @test GhpGhx_results["peak_combined_heatpump_thermal_ton"] ≈ 466.632 rtol=0.02    

end

@testset "WWHP system sizing" begin

    inputs_dict = JSON.parsefile("inputs/inputs_base.json")
    cooling_cop_map_df = CSV.read("inputs/wwhp_cooling_cop_map.csv", DataFrame)
    heating_cop_map_df = CSV.read("inputs/wwhp_heating_cop_map.csv", DataFrame)

    # Generate a "records" style dictionary from the 
    cooling_cop_map_list = []
    heating_cop_map_list = []
    for i in eachrow(cooling_cop_map_df)
        dict_record = Dict(names(i)[k]=>i[names(i)[k]] for k in 1:length(names(i)))
        push!(cooling_cop_map_list, dict_record)
    end
    for i in eachrow(heating_cop_map_df)
        dict_record = Dict(names(i)[k]=>i[names(i)[k]] for k in 1:length(names(i)))
        push!(heating_cop_map_list, dict_record)
    end
    inputs_dict["wwhp_cop_map_eft_cooling"] = cooling_cop_map_list
    inputs_dict["wwhp_cop_map_eft_heating"] = heating_cop_map_list

    inputs_dict["heat_pump_configuration"] = "WWHP"
    inputs_dict["wwhp_cooling_setpoint_f"] = 45
    inputs_dict["wwhp_heating_setpoint_f"] = 100

    # Heating and cooling loads
    inputs_dict["heating_thermal_load_mmbtu_per_hr"] = ones(8760) * 2
    inputs_dict["cooling_thermal_load_ton"] = ones(8760) * 300

    # TODO Add PVWatts API call for dry bulb temperature based on lat/long
    inputs_dict["ambient_temperature_f"] = ones(8760) * 59.0

    # TODO Add ground_k_by_climate_zone lookup functionality from REopt.jl and need lat/long too
    inputs_dict["ground_thermal_conductivity_btu_per_hr_ft_f"] = 1.2

    @info "Starting GhpGhx" #with timeout of $(timeout) seconds..."
    results, inputs_params = GhpGhx.ghp_model(inputs_dict)
    # Create a dictionary of the results data needed for REopt
    GhpGhx_results = GhpGhx.get_results_for_reopt(results, inputs_params)
    @info "GhpGhx model solved" #with status $(results["status"])."

    @test typeof(results) == GhpGhx.ResultsStruct
    @test GhpGhx_results["number_of_boreholes"] ≈ 1283 rtol=0.02
    @test GhpGhx_results["peak_cooling_heatpump_thermal_ton"] ≈ 300.0 rtol=0.02 
    @test GhpGhx_results["peak_heating_heatpump_thermal_ton"] ≈ 166.632 rtol=0.02 
end