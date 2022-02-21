#=
Advanced thermodynamics calculations. Once the sub-module Thermodynamics
of ClimateMachine.jl becomes a separate package, those can be contributed
there.
=#
#########################################################################################
# Constants
#########################################################################################
module Constants
const Lhvap = 2.5e6    # Latent heat of vaporization (J / kg)
const Lhsub = 2.834e6   # Latent heat of sublimation (J / kg)
const Lhfus = Lhsub - Lhvap  # Latent heat of fusion (J / kg)
const cp = 1004.0     # specific heat at constant pressure for dry air (J / kg / K)
const Rd = 287.0         # gas constant for dry air (J / kg / K)
const R_over_cp = Rd / cp
const Rv = 461.5       # gas constant for water vapor (J / kg / K)
const cpv = 1875.0     # specific heat at constant pressure for water vapor (J / kg / K)
const eps = Rd / Rv
const g = 9.8          # gravitational acceleration (m / s**2)
const rho_w = 1000.0    # density of water (kg / m**3)
const cw = 4181.3      # specific heat of liquid water (J / kg / K)
end

#########################################################################################
# EIS/ECTEI
#########################################################################################
"""
```julia
estimated_cloud_top_entrainment_index(
    Tsfc, T700, qsfc, q700;
    Psfc = 1000.0, RH = 0.8
) → ECTEI (Kelvin)
```
Calculate the ECTEI according to eq.(3) of [^Kawai2017] using [estimated_inversion_strength](@ref)
and the specific humidity at surface and 700 hPa:

```math
ECTEI = EIS - \\beta (L/c_p) (q_{sfc} - q_{700})
```
with ``\\beta, L, c_p`` constants.

ECTEI is an improvement of EIS as a predictor of low cloud cover.

[^Kawai2017]: Interpretation of factors controlling low cloud cover and low cloud feedback
              using a unified predictive index,
              Kawai et al 2017, [DOI: 10.1175/JCLI-D-16-0825.1](https://journals.ametsoc.org/view/journals/clim/30/22/jcli-d-16-0825.1.xml)
"""
function estimated_cloud_top_entrainment_index(
        Tsfc, T700, qsfc, q700;
        Psfc = 1000.0, RH = 0.8
    )
    EIS = estimated_inversion_strength(Tsfc, T700; Psfc, RH)
    β = 0.23 # page 2 of paper
    ECTEI = @. EIS - β*(Constants.Lhvap/Constants.cp)*(qsfc - q700)
end

# `estimated_inversion_strength` and `lifting_condensation_level` are from Python's `climlab`, see here:
# https://climlab.readthedocs.io/en/latest/_modules/climlab/utils/thermo.html#estimated_inversion_strength
# It has MIT license: https://climlab.readthedocs.io/en/latest/license.html

"""
    estimated_inversion_strength(Tsfc, T700; Psfc = 1000.0, RH = 0.8) → EIS
Compute the Estimated Inversion Strength or EIS, from eq.(4) of [^Wood2006],
given `Tsfc` the surface temperature and `T700` the air temperature at 700 hPa
(in Kelvin).

EIS is a normalized measure of lower tropospheric stability acccounting for
temperature-dependence of the moist adiabat.

[^Wood2006]: On the relationship between stratiform low cloud cover and lower-tropospheric stability,
             Wood and Bretheron, [DOI: 10.1175/JCLI3988.1](https://journals.ametsoc.org/view/journals/clim/19/24/jcli3988.1.xml)
"""
function estimated_inversion_strength(Tsfc, T700; Psfc = 1000.0, RH = 0.8)
    # Interpolate to 850 hPa
    T850 = @. (Tsfc+T700)/2.0;
    # Assume 80% relative humidity to compute LCL, appropriate for marine boundary layer
    LCL = lifting_condensation_level(Tsfc, RH)
    # Lower Tropospheric Stability (θ700 - θ0)
    LTS = lower_tropospheric_stability(Tsfc, T700, Psfc)
    #  Gammam  = -dθ/dz is the rate of potential temperature decrease along the moist adiabat
    #  in K / m
    Gammam = @. (Constants.g/Constants.cp*(1.0 - (1.0 + Constants.Lhvap*qsat(T850,850) / Constants.Rd / T850) /
             (1.0 + Constants.Lhvap^2 * qsat(T850, 850)/ Constants.cp/Constants.Rv/T850^2)))
    # Assume exponential decrease of pressure with scale height given by surface temperature
    z700 = @. (Constants.Rd*Tsfc/Constants.g)*log(Psfc/700.0)
    return @. LTS - Gammam*(z700 - LCL)
end

"""
    lower_tropospheric_stability(Tsfc, T700, Psfc = 1000.0)
By definition the difference of the potential temperature at 700hPa minus that at surface
pressure, by default = 1000.0 hPa.
"""
function lower_tropospheric_stability(Tsfc, T700, Psfc = 1000.0)
    LTS = potential_temperature(T700, 700) .- potential_temperature(Tsfc, Psfc)
end

"""
    lifting_condensation_level(T, RH) → LCL (meters)
Compute the Lifiting Condensation Level (LCL) for a given temperature and relative humidity.

This is height (relative to parcel height) at which the parcel would become saturated
during adiabatic ascent. It is based on an approximate formula from Bolton
(1980 MWR) as given by Romps (2017 JAS).
For an exact formula see Romps (2017 JAS), doi:10.1175/JAS-D-17-0102.1
"""
function lifting_condensation_level(T, RH)
    Tadj = @. T - 55.0  # in Kelvin
    return @. Constants.cp/Constants.g*(Tadj - (1/Tadj - log(RH)/2840.0)^(-1))
end

function potential_temperature(T, P)
    θ = @. T*(1000.0/P)^Constants.R_over_cp
    return θ
end

"""
    qsat(T,P) → qₛ (dimensionless)
Compute saturation specific humidity as function of temperature (in Kelvin)
and pressure (in hPa).
"""
function qsat(T,p)
    es = clausius_clapeyron(T)
    q = @. Constants.eps * es / (p - (1 - Constants.eps) * es )
    return q
end

"""
    clausius_clapeyron(T) → es (hPa)
Compute saturation vapor pressure as function of temperature T (in Kelvin).

Formula from Rogers and Yau "A Short Course in Cloud Physics" (Pergammon Press), p. 16
claimed to be accurate to within 0.1% between -30degC and 35 degC
Based on the paper by Bolton (1980, Monthly Weather Review).
"""
function clausius_clapeyron(T)
    Tcel = T .- 273.15
    es = @. 6.112 * exp(17.67*Tcel/(Tcel+243.5))
    return es
end
