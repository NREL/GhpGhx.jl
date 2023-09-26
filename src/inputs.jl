# Define some constant data structures for default inputs in InputsStruct
const ground_k_by_climate_zone = Dict([
    ("1A", 1.029),
    ("2A", 1.348),
    ("2B", 0.917),
    ("3A", 1.243),
    ("3B", 1.364),
    ("3C", 1.117),
    ("4A", 1.023),
    ("4B", 0.972),
    ("4C", 1.418),
    ("5A", 1.726),
    ("5B", 1.177),
    ("6A", 0.977),
    ("6B", 0.981),
    ("7", 1.271),
    ("8", 1.189)
])

const default_cop_map_list = [
    Dict{String, Real}("eft" => 20, "cool_cop" => 11.023, "heat_cop" => 3.351)
    Dict{String, Real}("eft" => 30, "cool_cop" => 11.023, "heat_cop" => 3.639)
    Dict{String, Real}("eft" => 40, "cool_cop" => 11.023, "heat_cop" => 4.161)
    Dict{String, Real}("eft" => 50, "cool_cop" => 10.481, "heat_cop" => 4.681)
    Dict{String, Real}("eft" => 60, "cool_cop" => 9.168, "heat_cop" => 5.081) 
    Dict{String, Real}("eft" => 70, "cool_cop" => 7.263, "heat_cop" => 5.678) 
    Dict{String, Real}("eft" => 80, "cool_cop" => 5.826, "heat_cop" => 6.047) 
    Dict{String, Real}("eft" => 90, "cool_cop" => 4.803, "heat_cop" => 6.341) 
    Dict{String, Real}("eft" => 100, "cool_cop" => 3.9, "heat_cop" => 6.341)  
    Dict{String, Real}("eft" => 110, "cool_cop" => 3.279, "heat_cop" => 6.341)
    Dict{String, Real}("eft" => 120, "cool_cop" => 2.707, "heat_cop" => 6.341)
]

const default_wwhp_cooling_cop_map_list = [
    Dict{String, Real}("eft" => 20, "40" => 11.023, "50" => 11.023, "60" => 11.023)
    Dict{String, Real}("eft" => 30, "40" => 11.023, "50" => 11.023, "60" => 11.023)
    Dict{String, Real}("eft" => 40, "40" => 11.023, "50" => 11.023, "60" => 11.023)
    Dict{String, Real}("eft" => 50, "40" => 10.481, "50" => 10.481, "60" => 10.481)
    Dict{String, Real}("eft" => 60, "40" => 9.168, "50" => 9.168, "60" => 9.1683) 
    Dict{String, Real}("eft" => 70, "40" => 7.263, "50" => 7.263, "60" => 7.263) 
    Dict{String, Real}("eft" => 80, "40" => 5.826, "50" => 5.826, "60" => 5.8263) 
    Dict{String, Real}("eft" => 90, "40" => 4.803, "50" => 4.803, "60" => 4.803) 
    Dict{String, Real}("eft" => 100, "40" => 3.9, "50" => 3.9, "60" => 3.9)  
    Dict{String, Real}("eft" => 110, "40" => 3.279, "50" => 3.279, "60" => 3.279)
    Dict{String, Real}("eft" => 120, "40" => 2.707, "50" => 2.707, "60" => 2.707)
]

const default_wwhp_heating_cop_map_list = [
    Dict{String, Real}("eft" => 20, "100" => 3.351, "110" => 3.351, "120" => 3.351)
    Dict{String, Real}("eft" => 30, "100" => 3.639, "110" => 3.639, "120" => 3.639)
    Dict{String, Real}("eft" => 40, "100" => 4.161, "110" => 4.161, "120" => 4.161)
    Dict{String, Real}("eft" => 50, "100" => 4.681, "110" => 4.681, "120" => 4.681)
    Dict{String, Real}("eft" => 60, "100" => 5.081, "110" => 5.081, "120" => 5.081) 
    Dict{String, Real}("eft" => 70, "100" => 5.678, "110" => 5.678, "120" => 5.678) 
    Dict{String, Real}("eft" => 80, "100" => 6.047, "110" => 6.047, "120" => 6.047) 
    Dict{String, Real}("eft" => 90, "100" => 6.341, "110" => 6.341, "120" => 6.341) 
    Dict{String, Real}("eft" => 100, "100" => 6.341, "110" => 6.341, "120" => 6.341)  
    Dict{String, Real}("eft" => 110, "100" => 6.341, "110" => 6.341, "120" => 6.341)
    Dict{String, Real}("eft" => 120, "100" => 6.341, "110" => 6.341, "120" => 6.341)
]

"""
    InputsStruct

This struct defines the inputs for the GhpGhx module
"""
Base.@kwdef mutable struct InputsStruct
    ##### These are the exact /GhpGhx POST names from the API #####
    # Parameters
    heat_pump_configuration::String  = "WSHP"  # "WSHP" or "WWHP"
    borehole_depth_ft::Float64 = 400.0
    ghx_header_depth_ft::Float64 = 4.0
    borehole_spacing_ft::Float64 = 20.0
    borehole_diameter_inch::Float64 = 5.0
    borehole_spacing_type::String  = "rectangular"  # "rectangular" or "hexagonal"
    ghx_pipe_outer_diameter_inch::Float64 = 1.66
    ghx_pipe_wall_thickness_inch::Float64 = 0.16
    ghx_pipe_thermal_conductivity_btu_per_hr_ft_f::Float64 = 0.25
    ghx_shank_space_inch::Float64 = 2.5
    ground_thermal_conductivity_btu_per_hr_ft_f::Float64 = NaN  # Default depends on climate zone
    ground_mass_density_lb_per_ft3::Float64 = 162.3
    ground_specific_heat_btu_per_lb_f::Float64 = 0.211
    grout_thermal_conductivity_btu_per_hr_ft_f::Float64 = 1.0
    ghx_fluid_specific_heat_btu_per_lb_f::Float64 = 1.0
    ghx_fluid_mass_density_lb_per_ft3::Float64 = 62.4
    ghx_fluid_thermal_conductivity_btu_per_hr_ft_f::Float64 = 0.34
    ghx_fluid_dynamic_viscosity_lbm_per_ft_hr::Float64 = 2.75
    ghx_fluid_flow_rate_gpm_per_ton::Float64 = 2.5
    ghx_pump_power_watt_per_gpm::Float64 = 15.0
    ghx_pump_min_speed_fraction::Float64 = 0.1
    ghx_pump_power_exponent::Float64 = 2.2
    max_eft_allowable_f::Float64 = 104.0
    min_eft_allowable_f::Float64 = 23.0
    hybrid_sizing_flag::Float64 = 1.0
    is_heating_electric::Bool = true
    aux_heater_thermal_efficiency::Float64 = 0.98
    aux_cooler_energy_use_intensity_kwe_per_kwt::Float64 = 0.02
    hybrid_ghx_sizing_method::String = "None"
    hybrid_ghx_sizing_fraction::Float64 = 0.6
    wwhp_cooling_setpoint_f::Float64 = 55.0
    wwhp_heating_setpoint_f::Float64 = 110.0

    # Centralized GHP
    wwhp_heating_pump_fluid_flow_rate_gpm_per_ton::Float64 = 3.0
    wwhp_cooling_pump_fluid_flow_rate_gpm_per_ton::Float64 = 3.0
    wwhp_heating_pump_power_watt_per_gpm::Float64 = 15.0
    wwhp_cooling_pump_power_watt_per_gpm::Float64 = 15.0
    wwhp_heating_pump_min_speed_fraction::Float64 = 0.1
    wwhp_cooling_pump_min_speed_fraction::Float64 = 0.1
    wwhp_heating_pump_power_exponent::Float64 = 2.2
    wwhp_cooling_pump_power_exponent::Float64 = 2.2

    # Array/Dict inputs
    heating_thermal_load_mmbtu_per_hr::Array{Float64,1} = Float64[]
    cooling_thermal_load_ton::Array{Float64,1} = Float64[]
    ambient_temperature_f::Array{Float64,1} = Float64[]
    cop_map_eft_heating_cooling::Array{Any,1} = Dict[]
    wwhp_cop_map_eft_heating::Array{Any,1} = Dict[]
    wwhp_cop_map_eft_cooling::Array{Any,1} = Dict[]
    
    # Model Settings
    hybrid_auto_ghx_sizing_flag::Bool = false # updates simulation_years and max_sizing_iterations for auto ghx sizing.
    simulation_years::Int64 = 25  # Number of years for GHP-GHX model
    solver_eft_tolerance_f::Float64 = 2.0  # Tolerance for the EFT error to accept a GHX sizing solution
    solver_eft_tolerance::Float64 = 2.0 / 1.8  # Convert to degC
    ghx_model::String = "TESS" # "TESS" or "DST"
    dst_ghx_timesteps_per_hour::Int64 = 12
    tess_ghx_minimum_timesteps_per_hour::Int64 = 1
    max_sizing_iterations::Int64 = 15
    init_sizing_factor_ft_per_peak_ton::Float64 = 75 #246.1
    
    ##### These are the variable names used in the GhpGhx, kept from TESS ######
    # TODO eventually just use the API names above in GhpGhx to remove redundancy (BUT WOULD HAVE TO DEAL WITH UNITS CONVERSION STILL)  
    I_Configuration::Int64 = 1  #!1=Decentralized WSHPs, 2=Decentralized WWHPs, 3=Centralized WWHPs, 4=Centralized WWHPs with HRC
    Depth_Bores::Float64 = NaN #!Depth of the boreholes (ft)
    Depth_Header::Float64 = NaN  #!Depth of the ground heat exchanger headers below grade (ft)
    BoreSpacing::Float64 = NaN  #!Spacing between boreholes (ft)
    SpacingType::Int64 = 1  # Spacing type (1=rectangular, 2 = hexagonal)
    Radius_Bores::Float64 = NaN  #!Radius of the boreholes (in)
    Ro_Pipe::Float64 = NaN  #!Outer radius of individual u-tube pipe (in)
    Ri_Pipe::Float64 = NaN  #!Inner radius of individual u-tube pipe (in)
    K_Pipe::Float64 = NaN  #!Thermal conductivity of pipe material (BTU/h.ft.F)
    Center_Center_Distance::Float64 = NaN  #!Distance between centers of upwards and downwards u-tube legs (in)
    K_Soil::Float64 = NaN  #!Thermal conductivity of the soil (BTU/h.ft.F)
    Rho_Soil::Float64 = NaN  #!Density of the soil (lbm/ft3)
    Cp_Soil::Float64 = NaN  #!Specific heat of the soil (BTU/lbm.F)
    K_Grout::Float64 = NaN  #!Thermal conductivity of the borehole grout (BTU/h.ft.F)
    T_Ground::Float64 = NaN  #!Average soil surface temperature over the year (F)
    Tamp_Ground::Float64 = NaN  #!Amplitude of soil surface temperature over the year (Delta_F)
    DayMin_Surface::Float64 = NaN  #!Day of minimum soil surface temperature (day)
    Cp_GHXFluid::Float64 = NaN #!Specific heat of the ground heat exchanger fluid (BTU/lbm.F)
    Rho_GHXFluid::Float64 = NaN  #!Density of the ground heat exchanger fluid (lbm/ft3)
    K_GHXFluid::Float64 = NaN  #!Thermal conductivity of the ground heat exchanger fluid (BTU/h.ft.F)
    Mu_GHXFluid::Float64 = NaN  #!Dynamic viscosity of the ground heat exchanger fluid (lbm/ft.h)
    GPMperTon_WSHP::Float64 = NaN  #!Nominal flow rate of the water source heat pumps (gpm/ton)
    WattPerGPM_GHXPump::Float64 = NaN  #!Nominal ground loop pump power (Watt/gpm)
    fMin_VSP_GHXPump::Float64 = NaN  #!Minimum ground heat exchanger pump speed for variable speed option (set to 1 if constant speed pumps)
    Exponent_GHXPump::Float64 = NaN  #!Exponent for relationship between ground heat exchanger pump power and ground heat exchanger pump flow rate
    Tmax_Sizing::Float64 = NaN  #!Maximum allowable return fluid temperature from the ground loop (F)
    Tmin_Sizing::Float64 = NaN  #!Minimum allowable return fluid temperature from the ground loop (F)    
    f_HybridSize::Float64 = NaN #Flag for hybrid GHX sizing (Modes - size for heating: -2; size for cooling: -1, size as a fraction of non-hybrid GHX: 0-2)

    # Additional parameters and site data
    HeatingThermalLoadKW::Array{Float64, 1} = Float64[]  # Heating thermal load to be served by GHP
    CoolingThermalLoadKW::Array{Float64, 1} = Float64[]  # Cooling thermal load to be served by GHP
    AmbientTemperature::Array{Float64, 1} = Float64[]  # Dry-bulb outdoor air temperature in degrees Fahrenheit 
    HeatPumpCOPMap::Matrix{Float64}  # Includes heating and cooling heat pump COP versus EWT (degF)
    HeatingHeatPumpCOPMap::Matrix{Float64}  # Heating heat pump COP versus EWT (degF)
    CoolingHeatPumpCOPMap::Matrix{Float64}  # Cooling heat pump COP versus EWT (degF)

    # These are defined based on the above inputs and processed in the function below
    TON_TO_KW::Float64 = NaN # [kw/ton]
    MMBTU_TO_KWH::Float64 = NaN  # [kwh/MMBtu]
    METER_TO_FEET::Float64 = NaN  # [ft/m]

    PeakTons_WSHP_H::Float64 = NaN
    PeakTons_WSHP_C::Float64 = NaN
    PeakTons_WSHP_GHX::Float64 = NaN
    X_init::Float64 = NaN
    N_Series::Int64 = 1
    N_Radial::Int64 = 1
    N_Vertical::Int64 = 1
    RhoCp_Soil::Float64 = NaN
    DayMin_DST::Float64 = NaN
    GPM_GHXPump::Float64 = NaN
    Prated_GHXPump::Float64 = NaN
    LPS_GHXPump::Float64 = NaN
    Mdot_GHXPump::Float64 = NaN

    # Centralized GHP - WWHP
    GPMperTon_WWHP_H::Float64 = NaN  #!Nominal flow rate of the water-water heating heat pumps from the GHX (gpm/ton) 
    GPMperTon_WWHP_C::Float64 = NaN  #!Nominal flow rate of the water-water cooling heat pumps from the GHX (gpm/ton)
    WattPerGPM_HeatingPumps::Float64 = NaN #!Nominal pump power for delivery of hot water to the loads (WWHP systems) (Watt/gpm)
    WattPerGPM_CoolingPumps::Float64 = NaN #!Nominal pump power for delivery of chilled water to the loads (WWHP systems) (Watt/gpm)
    fmin_VSP_HeatingPumps::Float64 = NaN  #!Minimum heating loop pump speed for variable speed option (set to 1 if constant speed pumps)
    fmin_VSP_CoolingPumps::Float64 = NaN  #!Minimum cooling pump speed for variable speed option (set to 1 if constant speed pumps)
    Exponent_HeatingPumps::Float64 = NaN  #!Exponent for relationship between heating loop pump power and heating loop pump flow rate
    Exponent_CoolingPumps::Float64 = NaN  #!Exponent for relationship between cooling loop pump power and cooling loop pump flow rate

    PeakTons_WWHP_H::Float64 = NaN
    PeakTons_WWHP_C::Float64 = NaN
    PeakTons_WWHP_GHX::Float64 = NaN
    fPeak_WWHP_H::Float64 = NaN

    GPM_HeatingPumps::Float64 = NaN
    GPM_CoolingPumps::Float64 = NaN
    Prated_HeatingPumps::Float64 = NaN
    Prated_CoolingPumps::Float64 = NaN
end


"""
    InputsProcess(d::Dict)

Performs unit conversions, name conversions, and additional processing of inputs.

"""
function InputsProcess(d::Dict)   
    # Convert all Dict key strings to symbols which is required for kwargs inputs of InputsStruct
    d = dictkeys_tosymbols(d)

    d[:HeatPumpCOPMap] = zeros(Float64, 1, 1)
    d[:HeatingHeatPumpCOPMap] = zeros(Float64, 1, 1)
    d[:CoolingHeatPumpCOPMap] = zeros(Float64, 1, 1)

    # Instantiate the mutable struct for assigning default values from Base.@kwdef and allows processing/modifying
    d = InputsStruct(; d...)

    if d.hybrid_auto_ghx_sizing_flag
        @info "Running GhpGhx for automatic hybrid GHX sizing. Simulation years \
         is 2 and maximum sizing iterations is 1"
        d.simulation_years = 2
        d.max_sizing_iterations = 1
    end

    # Heat pump configuration
    if d.heat_pump_configuration == "WSHP"
        d.I_Configuration = 1
    elseif d.heat_pump_configuration == "WWHP"
        d.I_Configuration = 3
    else
        print("Unknown heat pump configuration entered.")
    end
    
    if d.I_Configuration == 1
        # Load in default COP map, if not input, which is a NON-keyword argument, so required for InputsStruct instantiation
        if isempty(d.cop_map_eft_heating_cooling)
            cop_map_list = deepcopy(default_cop_map_list)
        else
            cop_map_list = d.cop_map_eft_heating_cooling
        end
        # Convert COP map list_of_dict to Matrix{Float64} (alias for Array{Float64, 2})
        d.HeatPumpCOPMap = zeros(Float64, length(cop_map_list), 3)
        for i in eachindex(cop_map_list)
            d.HeatPumpCOPMap[i,1] = cop_map_list[i]["eft"]
            d.HeatPumpCOPMap[i,2] = cop_map_list[i]["heat_cop"]
            d.HeatPumpCOPMap[i,3] = cop_map_list[i]["cool_cop"]
        end  

    elseif d.I_Configuration == 3
        
        if isempty(d.wwhp_cop_map_eft_heating)
            heating_cop_map_list = deepcopy(default_wwhp_heating_cop_map_list)
        else
            heating_cop_map_list = d.wwhp_cop_map_eft_heating
        end
        if isempty(d.wwhp_cop_map_eft_cooling)
            cooling_cop_map_list = deepcopy(default_wwhp_cooling_cop_map_list)
        else
            cooling_cop_map_list = d.wwhp_cop_map_eft_cooling
        end

        heating_COPs, heating_EFTs = get_wwhp_cop_matrix(heating_cop_map_list, d.wwhp_heating_setpoint_f)
        cooling_COPs, cooling_EFTs = get_wwhp_cop_matrix(cooling_cop_map_list, d.wwhp_cooling_setpoint_f)

        d.HeatingHeatPumpCOPMap = zeros(Float64, length(heating_cop_map_list), 2)
        d.CoolingHeatPumpCOPMap = zeros(Float64, length(cooling_cop_map_list), 2)

        d.HeatingHeatPumpCOPMap[:, 1] = heating_EFTs
        d.HeatingHeatPumpCOPMap[:, 2] = heating_COPs
        d.CoolingHeatPumpCOPMap[:, 1] = cooling_EFTs
        d.CoolingHeatPumpCOPMap[:, 2] = cooling_COPs
    end

    # Constants
    d.TON_TO_KW = 3.51685  # [kw/ton]
    d.MMBTU_TO_KWH = 293.0107  # [kwh/MMBtu]
    d.METER_TO_FEET= 3.28084  # [ft/m]

    # Hybrid Flag
    d.f_HybridSize = d.hybrid_sizing_flag

    # Convert API inputs to GhpGhx variable names, and units from English to SI
    d.Depth_Bores=  d.borehole_depth_ft/ d.METER_TO_FEET # [m]
    d.Depth_Header= d.ghx_header_depth_ft/ d.METER_TO_FEET # [m]
    d.BoreSpacing= d.borehole_spacing_ft/ d.METER_TO_FEET # [m]
    if d.borehole_spacing_type== "rectangular"
        d.SpacingType= 1
    else
        d.SpacingType= 2
    end
    d.Radius_Bores= d.borehole_diameter_inch/ 2.0 / 12.0 / d.METER_TO_FEET # [m]
    d.Ro_Pipe= d.ghx_pipe_outer_diameter_inch/ 2.0 / 12.0 / d.METER_TO_FEET # [m]
    d.Ri_Pipe= (d.ghx_pipe_outer_diameter_inch/ 2.0 - d.ghx_pipe_wall_thickness_inch) / 12.0 / d.METER_TO_FEET # [m]
    d.K_Pipe= d.ghx_pipe_thermal_conductivity_btu_per_hr_ft_f* 1.055 * d.METER_TO_FEET* 1.8  # [kJ/h.m.K]
    d.Center_Center_Distance= d.ghx_shank_space_inch/ 12.0 / d.METER_TO_FEET # [m]
    d.K_Soil= d.ground_thermal_conductivity_btu_per_hr_ft_f* 1.055 * d.METER_TO_FEET* 1.8  # [kJ/h.m]
    d.Rho_Soil= d.ground_mass_density_lb_per_ft3* 35.31467 / 2.20462  # [kg/m3]
    d.Cp_Soil= d.ground_specific_heat_btu_per_lb_f* 1.8 * 2.20462 * 1.055  # [kJ/kg.K]
    d.K_Grout= d.grout_thermal_conductivity_btu_per_hr_ft_f* 1.055 * d.METER_TO_FEET* 1.8  # [kJ/h.m.K   
    # See processing of ambient temperature for temperature parameters, below
    d.Cp_GHXFluid= d.ghx_fluid_specific_heat_btu_per_lb_f* 1.8 * 2.20462 * 1.055  # [kJ/kg.K]
    d.Rho_GHXFluid= d.ghx_fluid_mass_density_lb_per_ft3* 35.31467 / 2.20462  # [kg/m3]
    d.K_GHXFluid= d.ghx_fluid_thermal_conductivity_btu_per_hr_ft_f* 1.055 * d.METER_TO_FEET* 1.8  # [kJ/h.m.K]
    d.Mu_GHXFluid= d.ghx_fluid_dynamic_viscosity_lbm_per_ft_hr* d.METER_TO_FEET/ 2.20462  # [kg/h.m]
    d.GPMperTon_WSHP= d.ghx_fluid_flow_rate_gpm_per_ton
    d.WattPerGPM_GHXPump= d.ghx_pump_power_watt_per_gpm
    d.fMin_VSP_GHXPump= d.ghx_pump_min_speed_fraction
    d.Exponent_GHXPump= d.ghx_pump_power_exponent
    d.Tmax_Sizing= (d.max_eft_allowable_f- 32.0) / 1.8  # [C]
    d.Tmin_Sizing= (d.min_eft_allowable_f- 32.0) / 1.8  # [C]
    d.solver_eft_tolerance= d.solver_eft_tolerance_f / 1.8  # [C]

    # Convert array input units to SI
    d.HeatingThermalLoadKW= d.heating_thermal_load_mmbtu_per_hr* d.MMBTU_TO_KWH
    d.CoolingThermalLoadKW= d.cooling_thermal_load_ton* d.TON_TO_KW
    d.AmbientTemperature= (d.ambient_temperature_f.- 32.0) / 1.8  # [C]

    # Use AmbientTemperature to calculate other temperature parameters
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    avg_temp_month = zeros(12)
    for (mo,days) in enumerate(days_in_month)
        hour_start = sum(days_in_month[1:mo-1]) * 24 + 1  # mo=1 results in 0 + 1 = 1
        hour_end = sum(days_in_month[1:mo]) * 24
        avg_temp_month[mo] = sum(d.AmbientTemperature[hour_start:hour_end]) / (days * 24)
    end
    
    d.T_Ground= sum(d.AmbientTemperature) / (365 * 24)  # [C]
    d.Tamp_Ground= (maximum(avg_temp_month) - minimum(avg_temp_month)) / 2  # [C]
    d.DayMin_Surface= convert(Int64, round(argmin(d.AmbientTemperature) / 24))  # day of year

    # Find peak heating, cooling, and combined for initial sizing guess
    d.PeakTons_WSHP_H= maximum(d.HeatingThermalLoadKW) / d.TON_TO_KW
    d.PeakTons_WSHP_C= maximum(d.cooling_thermal_load_ton)
    d.PeakTons_WSHP_GHX= maximum(d.HeatingThermalLoadKW + d.CoolingThermalLoadKW) / d.TON_TO_KW
    d.X_init= d.init_sizing_factor_ft_per_peak_ton/ d.METER_TO_FEET* d.PeakTons_WSHP_GHX # [m]
    
    # Centralized GHP - WWHP 
    d.GPMperTon_WWHP_H = d.wwhp_heating_pump_fluid_flow_rate_gpm_per_ton
    d.GPMperTon_WWHP_C = d.wwhp_cooling_pump_fluid_flow_rate_gpm_per_ton
    d.WattPerGPM_HeatingPumps = d.wwhp_heating_pump_power_watt_per_gpm
    d.WattPerGPM_CoolingPumps = d.wwhp_cooling_pump_power_watt_per_gpm
    d.fmin_VSP_HeatingPumps = d.wwhp_heating_pump_min_speed_fraction
    d.fmin_VSP_CoolingPumps = d.wwhp_cooling_pump_min_speed_fraction
    d.Exponent_HeatingPumps = d.wwhp_heating_pump_power_exponent
    d.Exponent_CoolingPumps = d.wwhp_cooling_pump_power_exponent

    d.PeakTons_WWHP_H = maximum(d.HeatingThermalLoadKW) / d.TON_TO_KW
    d.PeakTons_WWHP_C = maximum(d.cooling_thermal_load_ton)
    d.PeakTons_WWHP_GHX = maximum(d.HeatingThermalLoadKW + d.CoolingThermalLoadKW) / d.TON_TO_KW
    idx = findmax(d.HeatingThermalLoadKW + d.CoolingThermalLoadKW)[2]
    d.fPeak_WWHP_H = (d.HeatingThermalLoadKW / d.TON_TO_KW)[idx] / (max(0.0001, (d.HeatingThermalLoadKW / d.TON_TO_KW)[idx] + d.cooling_thermal_load_ton[idx]))
    
    # Set some intermediate conditions
    d.N_Series= 1
    d.N_Radial= d.N_Series
    d.N_Vertical= 50
    d.RhoCp_Soil= d.Rho_Soil* d.Cp_Soil
    d.DayMin_DST= 270.0 - d.DayMin_Surface

    if d.I_Configuration == 1
        d.GPM_GHXPump = d.GPMperTon_WSHP* d.PeakTons_WSHP_GHX
    elseif d.I_Configuration == 3
        d.GPM_GHXPump = d.GPMperTon_WWHP_H * d.PeakTons_WWHP_GHX * d.fPeak_WWHP_H + d.GPMperTon_WWHP_C * d.PeakTons_WWHP_GHX * (1 - d.fPeak_WWHP_H)

        d.GPM_HeatingPumps = d.GPMperTon_WWHP_H * d.PeakTons_WWHP_H 
        d.GPM_CoolingPumps = d.GPMperTon_WWHP_C * d.PeakTons_WWHP_C 
        d.Prated_HeatingPumps = d.WattPerGPM_HeatingPumps * 3.6 * d.GPM_HeatingPumps
        d.Prated_CoolingPumps = d.WattPerGPM_CoolingPumps * 3.6 * d.GPM_CoolingPumps
    end

    d.Prated_GHXPump= d.WattPerGPM_GHXPump* 3.6 * d.GPM_GHXPump
    d.LPS_GHXPump= d.GPM_GHXPump/ 60 / 264.172 * 1000.0
    d.Mdot_GHXPump= d.GPM_GHXPump* 60 / 264.172 * d.Rho_GHXFluid

    # Return processed (modified) InputsStruct d
    return d 
end

function get_wwhp_cop_matrix(cop_map::Array, setpoint::Float64)
    
    cop_map_keys = deepcopy(cop_map[1])
    pop!(cop_map_keys, "eft")
    
    temperature_setpoints = sort(parse.(Int64, collect(keys(cop_map_keys))))
    
    HeatPumpCOPMap = zeros(Float64, length(cop_map), length(temperature_setpoints)+1)
    
    # Create heat pump COP matrix
    for i in eachindex(cop_map)
        HeatPumpCOPMap[i,1] = cop_map[i]["eft"]
        for (idx, tmp) in enumerate(temperature_setpoints)
            HeatPumpCOPMap[i, idx+1] = cop_map[i][string(tmp)]
        end 
    end 
    
    if setpoint <= temperature_setpoints[1]
        cop = [:, 2]
    elseif setpoint > temperature_setpoints[end]
        cop = HeatPumpCOPMap[:, end]     
    elseif setpoint in temperature_setpoints
        idx = findfirst(x->x == setpoint, temperature_setpoints)
        cop = HeatPumpCOPMap[:, idx+1]  
    else
        # Loop over all temperature points in HeatPumpCOPMap, but break out of loop if it is found
        for (index, tmp) in enumerate(temperature_setpoints[2:end])  # Omit first and last temp checks, done above    
            if setpoint > temperature_setpoints[index, 1] && setpoint <= tmp   
                cop = []
                for i in eachindex(HeatPumpCOPMap[:,1])
                    slope = (HeatPumpCOPMap[i, index+2] - HeatPumpCOPMap[i, index+1]) / (tmp - temperature_setpoints[index])
                    append!(cop, HeatPumpCOPMap[i, index+1] + (tmp - setpoint) * slope)
                end
                break
            end
        end
    end

    eft = HeatPumpCOPMap[:,1]
    
    return cop, eft
end


"""
    dictkeys_tosymbols(d::Dict)

Assists in using a Dict as an input to instantiate the struct.

"""
function dictkeys_tosymbols(d::Dict)
    d2 = Dict()
    for (k, v) in d
        d2[Symbol(k)] = v
    end
    return d2
end