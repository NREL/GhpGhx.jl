"""
    ResultsStruct

This struct gets updated in the GhpGhx.jl module and is used to record the results/outputs.
"""
Base.@kwdef mutable struct ResultsStruct
    # These get modified by the function below
    # Some of these values are required inputs, those that do not have "= value"
    X_Now::Array{Float64, 1}
    FX_Now::Array{Float64, 1}
    N_Bores::Array{Int64, 1}
    N_Bores_Final::Int64 = 0
    FX_Final::Float64 = 0.0
    Length_Boreholes::Float64 = 0.0
    StorageVolume::Float64 = 0.0
    Depth_DST::Float64 = 0.0
    Mdot_Borehole::Float64 = 0.0
    Radius_Field::Float64 = 0.0
    N_Circuits::Int64 = 0
    L_Circuit::Float64 = 0.0
    Mass_Circuit::Float64 = 0.0

    # These get reset each sizing iteration
    Tmax_GHX::Float64 = -9.99e20
    Tmin_GHX::Float64 = 9.99e20
    Power_GHXPump::Float64 = 0.0
    Qfluid_GHXPump::Float64 = 0.0
    Power_WSHP_H::Float64 = 0.0
    Power_WSHP_C::Float64 = 0.0
    LoadMet_WSHP_H::Float64 = 0.0
    LoadMet_WSHP_C::Float64 = 0.0
    Q_Rejected_Total::Float64 = 0.0
    Q_Absorbed_Total::Float64 = 0.0
    Q_GHX_Net::Float64 = 0.0
    Q_GHX_In::Float64 = 0.0
    Q_GHX_Top::Float64 = 0.0
    Q_GHX_Sides::Float64 = 0.0
    Q_GHX_Bottom::Float64 = 0.0
    Q_GHX_Stored::Float64 = 0.0
    Qf_GHX_Stored::Float64 = 0.0
    Q_AuxHeat_Total::Float64 = 0.0
    Q_AuxCool_Total::Float64 = 0.0
    Qf_GHXPump::Float64 = 0.0
    Qnet_GHX::Float64 = 0.0
    Qnet_HeatPumps::Float64 = 0.0
    Qnet_RunningTotal::Float64 = 0.0
    
    # Centralized GHP
    Power_HeatingPumps::Float64 = 0.0
    Qfluid_HeatingPumps::Float64 = 0.0
    Power_CoolingPumps::Float64 = 0.0
    Qfluid_CoolingPumps::Float64 = 0.0
    
    # Array results
    total_hours::Int64  # Required argument to be passed by InputsStruct
    P_GHXPump_Hourly::Array{Float64, 1} = zeros(total_hours)
    P_HeatingPumps_Hourly::Array{Float64, 1} = zeros(total_hours)
    P_CoolingPumps_Hourly::Array{Float64, 1} = zeros(total_hours)
    P_WSHPh_Hourly::Array{Float64, 1} = zeros(total_hours)
    P_WSHPc_Hourly::Array{Float64, 1} = zeros(total_hours)
    Qh_Hourly::Array{Float64, 1} = zeros(total_hours)
    Qc_Hourly::Array{Float64, 1} = zeros(total_hours)
    QauxHt_Hourly::Array{Float64, 1} = zeros(total_hours)
    QauxCl_Hourly::Array{Float64, 1} = zeros(total_hours)
    EWT::Array{Float64, 1} = zeros(total_hours)
    LWT::Array{Float64, 1} = zeros(total_hours)
    N_iterations::Int64 = 0
end

"""
    init_sizing!(r::ResultsStruct, p::InputsStruct, size_iter::Int64)

Performs the initial sizing for the GhpGhx.jl `size_borefield(p)` function.
"""
function init_sizing!(r::ResultsStruct, p::InputsStruct, size_iter::Int64)
    #Changes each sizing iteration
    r.N_Bores[size_iter] = floor(r.X_Now[size_iter] / p.Depth_Bores) + 1
    r.Length_Boreholes = max(r.X_Now[size_iter] / r.N_Bores[size_iter], 0.001)
    if p.SpacingType == 1
        r.StorageVolume = r.Length_Boreholes * r.N_Bores[size_iter] * pi * (p.BoreSpacing * 0.525)^2
    else
        r.StorageVolume = p.BoreSpacing^2 * r.Length_Boreholes * r.N_Bores[size_iter]
    end
    r.Depth_DST = 1.5 * r.Length_Boreholes
    r.Mdot_Borehole = p.Mdot_GHXPump / r.N_Bores[size_iter] * p.N_Series
    r.Radius_Field = (r.StorageVolume / (pi * r.Length_Boreholes))^0.5
    r.N_Circuits = r.N_Bores[size_iter] / p.N_Series
    r.L_Circuit = 2.0 * p.N_Series * r.Length_Boreholes + (p.N_Series - 1) / p.N_Series * r.Radius_Field
    r.Mass_Circuit = r.L_Circuit * pi * p.Ri_Pipe^2 * p.Rho_GHXFluid
    nothing
end

"""
    get_GhpGhx_results_for_reopt(r::ResultsStruct, p::InputsStruct)

Extract the required results for input to REopt /job: average hourly dispatch performance and other summary performance metrics.
"""
function get_results_for_reopt(r::ResultsStruct, p::InputsStruct)
    results_dict = Dict{Any,Any}()
    results_dict["number_of_boreholes"] = r.N_Bores_Final
    results_dict["length_boreholes_ft"] = round(r.Length_Boreholes * p.METER_TO_FEET, digits=1)
    results_dict["yearly_heating_heatpump_electric_consumption_series_kw"] = zeros(8760)
    results_dict["yearly_cooling_heatpump_electric_consumption_series_kw"] = zeros(8760)
    results_dict["yearly_ghx_pump_electric_consumption_series_kw"] = zeros(8760)
    results_dict["yearly_heat_pump_eft_series_f"] = zeros(8760)
    # Hybrid
    results_dict["yearly_aux_heater_thermal_production_series_mmbtu_per_hour"] = zeros(8760)
    yearly_aux_cooler_thermal_production_series_ton = zeros(8760)
    yearly_aux_heater_thermal_production_series_kwt = zeros(8760)
    results_dict["yearly_aux_cooler_thermal_production_series_kwt"] = zeros(8760)
    results_dict["yearly_aux_heater_electric_consumption_series_kw"] = zeros(8760)
    results_dict["yearly_aux_cooler_electric_consumption_series_kw"] = zeros(8760)
    results_dict["yearly_ghx_lft_series_f"] = zeros(8760)
    results_dict["end_of_year_ghx_lft_f"] = zeros(0)
    results_dict["max_yearly_ghx_lft_f"] = zeros(0)
    results_dict["min_yearly_ghx_lft_f"] = zeros(0)
    
    # Get average electric consumption
    for yr in 1:p.simulation_years
        results_dict["yearly_heating_heatpump_electric_consumption_series_kw"] += round.(r.P_WSHPh_Hourly[(yr-1)*8760+1:yr*8760] / p.simulation_years, digits=3)
        results_dict["yearly_cooling_heatpump_electric_consumption_series_kw"] += round.(r.P_WSHPc_Hourly[(yr-1)*8760+1:yr*8760] / p.simulation_years, digits=3)
        results_dict["yearly_ghx_pump_electric_consumption_series_kw"] += round.(r.P_GHXPump_Hourly[(yr-1)*8760+1:yr*8760] / p.simulation_years, digits=3)
        results_dict["yearly_heat_pump_eft_series_f"] += round.((r.EWT[(yr-1)*8760+1:yr*8760] * 1.8 .+ 32.0) / p.simulation_years, digits=3)
    end

    # Hybrid - get average auxiliary heater/cooler production
    for yr in 1:p.simulation_years
        results_dict["yearly_ghx_lft_series_f"] += round.((r.LWT[(yr-1)*8760+1:yr*8760] * 1.8 .+ 32.0) / p.simulation_years, digits=3)
        yearly_aux_heater_thermal_production_series_kwt += round.(r.QauxHt_Hourly[(yr-1)*8760+1:yr*8760] / p.simulation_years, digits=3)
        results_dict["yearly_aux_cooler_thermal_production_series_kwt"] += round.(r.QauxCl_Hourly[(yr-1)*8760+1:yr*8760] / p.simulation_years, digits=3)
        append!(results_dict["end_of_year_ghx_lft_f"], round.((r.LWT[yr*8760] * 1.8 .+ 32.0), digits=3))
        append!(results_dict["max_yearly_ghx_lft_f"], round.(maximum(r.LWT[(yr-1)*8760+1:yr*8760]* 1.8 .+ 32.0), digits=3))
        append!(results_dict["min_yearly_ghx_lft_f"], round.(minimum(r.LWT[(yr-1)*8760+1:yr*8760] * 1.8 .+ 32.0), digits=3))
    end

    results_dict["yearly_aux_heater_thermal_production_series_mmbtu_per_hour"] = yearly_aux_heater_thermal_production_series_kwt / p.MMBTU_TO_KWH
    yearly_aux_cooler_thermal_production_series_ton = results_dict["yearly_aux_cooler_thermal_production_series_kwt"]  / p.TON_TO_KW
    results_dict["peak_aux_heater_thermal_production_mmbtu_per_hour"] = round(maximum(results_dict["yearly_aux_heater_thermal_production_series_mmbtu_per_hour"]), digits=3)
    results_dict["peak_aux_cooler_thermal_production_ton"] = round(maximum(yearly_aux_cooler_thermal_production_series_ton), digits=3)

    # Elec consumption of auxiliary cooler and/or heater
    results_dict["yearly_aux_cooler_electric_consumption_series_kw"] = results_dict["yearly_aux_cooler_thermal_production_series_kwt"] * p.aux_cooler_energy_use_intensity_kwe_per_kwt
    if p.is_heating_electric
        results_dict["yearly_aux_heater_electric_consumption_series_kw"] = yearly_aux_heater_thermal_production_series_kwt / p.aux_heater_thermal_efficiency
    end
    
    results_dict["annual_aux_heater_electric_consumption_kwh"] = round(sum(results_dict["yearly_aux_heater_electric_consumption_series_kw"]), digits=3) 
    results_dict["annual_aux_cooler_electric_consumption_kwh"] = round(sum(results_dict["yearly_aux_cooler_electric_consumption_series_kw"]), digits=3)

    results_dict["aux_heat_exchange_unit_type"] = "None"
    if results_dict["annual_aux_heater_electric_consumption_kwh"] > 0.1
        results_dict["aux_heat_exchange_unit_type"] = "Heater"
    elseif results_dict["annual_aux_cooler_electric_consumption_kwh"] > 0.1
        results_dict["aux_heat_exchange_unit_type"] = "Cooler"
    end

    results_dict["ghx_soln_number_of_iterations"] = r.N_iterations

    results_dict["yearly_total_electric_consumption_series_kw"] = 
        results_dict["yearly_heating_heatpump_electric_consumption_series_kw"] + 
        results_dict["yearly_cooling_heatpump_electric_consumption_series_kw"] + 
        results_dict["yearly_ghx_pump_electric_consumption_series_kw"] + 
        results_dict["yearly_aux_heater_electric_consumption_series_kw"] + 
        results_dict["yearly_aux_cooler_electric_consumption_series_kw"]
    results_dict["yearly_total_electric_consumption_kwh"] = round(sum(results_dict["yearly_total_electric_consumption_series_kw"]), digits=1)
    results_dict["peak_heating_heatpump_thermal_ton"] = round(p.PeakTons_WSHP_H, digits=3)
    results_dict["peak_cooling_heatpump_thermal_ton"] = round(p.PeakTons_WSHP_C, digits=3)
    results_dict["peak_combined_heatpump_thermal_ton"] = round(p.PeakTons_WSHP_GHX, digits=3)
    results_dict["max_eft_f"] = round(maximum(r.EWT) * 1.8 + 32.0)
    results_dict["min_eft_f"] = round(minimum(r.EWT) * 1.8 + 32.0)
    # Calculate average COP for heating and cooling; estimate allocation of pump power to heating and cooling by thermal energy served
    heating_thermal_kwh = sum(p.HeatingThermalLoadKW)
    cooling_thermal_kwh = sum(p.CoolingThermalLoadKW)
    heating_thermal_frac = heating_thermal_kwh / (heating_thermal_kwh + cooling_thermal_kwh)
    cooling_thermal_frac = cooling_thermal_kwh / (heating_thermal_kwh + cooling_thermal_kwh)
    ghx_pump_electric_kwh = sum(results_dict["yearly_ghx_pump_electric_consumption_series_kw"])
    heating_ghx_pump_electric_kwh = heating_thermal_frac * ghx_pump_electric_kwh
    cooling_ghx_pump_electric_kwh = cooling_thermal_frac * ghx_pump_electric_kwh
    heating_heatpump_electric_kwh = sum(results_dict["yearly_heating_heatpump_electric_consumption_series_kw"])
    cooling_heatpump_electric_kwh = sum(results_dict["yearly_cooling_heatpump_electric_consumption_series_kw"])
    results_dict["heating_cop_avg"] = round(heating_thermal_kwh / (heating_heatpump_electric_kwh + heating_ghx_pump_electric_kwh), digits=3)
    results_dict["cooling_cop_avg"] = round(cooling_thermal_kwh / (cooling_heatpump_electric_kwh + cooling_ghx_pump_electric_kwh), digits=3)
    results_dict["solved_eft_error_f"] = round(r.FX_Final * 1.8, digits=3)
    return results_dict
end
