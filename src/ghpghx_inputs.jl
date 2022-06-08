# Warning: I've seen notes of Base.@kwdef not being officially supported (not documented)
"""
    InputsStruct

This struct defines the inputs for the GhpGhx module
"""
Base.@kwdef mutable struct InputsStruct
    ##### These are the exact /GhpGhx POST names from the API #####
    # Parameters
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
    
    # Array/Dict inputs
    heating_thermal_load_mmbtu_per_hr::Array{Float64,1} = Float64[]
    cooling_thermal_load_ton::Array{Float64,1} = Float64[]
    ambient_temperature_f::Array{Float64,1} = Float64[]
    cop_map_eft_heating_cooling::Array{Any,1} = Dict[]

    # Model Settings
    simulation_years::Int64 = 25  # Number of years for GHP-GHX model
    solver_eft_tolerance_f::Float64 = 2.0  # Tolerance for the EFT error to accept a GHX sizing solution
    solver_eft_tolerance::Float64 = 2.0 / 1.8  # Convert to degC
    ghx_model::String = "TESS" # "TESS" or "DST"
    dst_ghx_timesteps_per_hour::Int64 = 12
    tess_ghx_minimum_timesteps_per_hour::Int64 = 1
    max_sizing_iterations::Int64 = 15
    init_sizing_factor_ft_per_peak_ton::Float64 = 246.1
    
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

    # Additional parameters and site data
    HeatingThermalLoadKW::Array{Float64, 1} = Float64[]  # Heating thermal load to be served by GHP
    CoolingThermalLoadKW::Array{Float64, 1} = Float64[]  # Cooling thermal load to be served by GHP
    AmbientTemperature::Array{Float64, 1} = Float64[]  # Dry-bulb outdoor air temperature in degrees Fahrenheit 
    HeatPumpCOPMap::Matrix{Float64}  # Includes heating and cooling heat pump COP versus EWT (degF)

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
end

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

"""
    InputsProcess(d::Dict)

Performs unit conversions, name conversions, and additional processing of inputs.

"""
function InputsProcess(d::Dict)   
    # Convert all Dict key strings to symbols which is required for kwargs inputs of InputsStruct
    d = dictkeys_tosymbols(d)

    # Load in default COP map, if not input, which is a NON-keyword argument, so required for InputsStruct instantiation
    if !haskey(d, :cop_map_eft_heating_cooling) || isempty(d[:cop_map_eft_heating_cooling])
        cop_map_df = CSV.read("test/inputs/cop_map.csv", DataFrame)
        # Generate a "records" style dictionary from the 
        cop_map_list = []
        for i in eachrow(cop_map_df)
            dict_record = Dict(names(i)[k]=>i[names(i)[k]] for k in 1:length(names(i)))
            push!(cop_map_list, dict_record)
        end
    else
        cop_map_list = d[:cop_map_eft_heating_cooling]
    end
    # Convert COP map list_of_dict to Matrix{Float64} (alias for Array{Float64, 2})
    d[:HeatPumpCOPMap] = zeros(Float64, length(cop_map_list), 3)
    for i in 1:length(d[:cop_map_eft_heating_cooling])
        d[:HeatPumpCOPMap][i,1] = cop_map_list[i]["eft"]
        d[:HeatPumpCOPMap][i,2] = cop_map_list[i]["heat_cop"]
        d[:HeatPumpCOPMap][i,3] = cop_map_list[i]["cool_cop"]
    end    
    
    # Instantiate the mutable struct for assigning default values from Base.@kwdef and allows processing/modifying
    d = InputsStruct(; d...)

    # Constants
    d.TON_TO_KW = 3.5169  # [kw/ton]
    d.MMBTU_TO_KWH = 293.07  # [kwh/MMBtu]
    d.METER_TO_FEET= 3.28084  # [ft/m]

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
    d.Tamp_Ground=  (maximum(avg_temp_month) - minimum(avg_temp_month)) / 2  # [C]
    d.DayMin_Surface= convert(Int64, round(argmin(d.AmbientTemperature) / 24))  # day of year

    # Find peak heating, cooling, and combined for initial sizing guess
    d.PeakTons_WSHP_H= maximum(d.HeatingThermalLoadKW) / d.TON_TO_KW
    d.PeakTons_WSHP_C= maximum(d.cooling_thermal_load_ton)
    d.PeakTons_WSHP_GHX= maximum(d.HeatingThermalLoadKW+ d.CoolingThermalLoadKW) / d.TON_TO_KW
    d.X_init= d.init_sizing_factor_ft_per_peak_ton/ d.METER_TO_FEET* d.PeakTons_WSHP_GHX # [m]
    
    # Set some intermediate conditions
    d.N_Series= 1
    d.N_Radial= d.N_Series
    d.N_Vertical= 50
    d.RhoCp_Soil= d.Rho_Soil* d.Cp_Soil
    d.DayMin_DST= 270.0 - d.DayMin_Surface
    d.GPM_GHXPump= d.GPMperTon_WSHP* d.PeakTons_WSHP_GHX
    d.Prated_GHXPump= d.WattPerGPM_GHXPump* 3.6 * d.GPM_GHXPump
    d.LPS_GHXPump= d.GPM_GHXPump/ 60 / 264.172 * 1000.0
    d.Mdot_GHXPump= d.GPM_GHXPump* 60 / 264.172 * d.Rho_GHXFluid
    
    # Return processed (modified) InputsStruct d
    return d 
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