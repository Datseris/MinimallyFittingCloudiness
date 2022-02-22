# In this file we provide functions that create definitions of
# cloudiness and return them as spatiotemporal data
using ClimateBase

"According to paper by Datseris & Stevens."
function effective_cloud_albedo(EBAF)

    R = ncread(EBAF, "toa_sw_all_mon")
    K = ncread(EBAF, "toa_sw_clr_t_mon")
    I = ncread(EBAF, "solar_mon")
    F = ncread(EBAF, "cldarea_total_daynight_mon") ./ 100
    T = ncread(EBAF, "cldtau_total_day_mon")
    Tcorrected = sinusoidal_continuation(T, [1, 2.0]; Tmin = 0)
    C1 = _cloud_effective_albedo(F, Tcorrected, Float32(0.9))
    C1_timeavg = timemean(C1)

    # use `\:arrow_down:` and tab for arrows
    F_s_⬆ = ncread(EBAF, "sfc_sw_up_all_mon")
    F_s_⬇ = ncread(EBAF, "sfc_sw_down_all_mon")
    F_s_⬆_K = ncread(EBAF, "sfc_sw_up_clr_t_mon")
    F_s_⬇_K = ncread(EBAF, "sfc_sw_down_clr_t_mon");

    l = size(F_s_⬆, Time)
    argsall = timemean.((I[Time(1:l)], R[Time(1:l)], F_s_⬆, F_s_⬇))
    argsclr = timemean.((I[Time(1:l)], K[Time(1:l)], F_s_⬆_K, F_s_⬇_K))
    α_ATM, α_SFC = surface_atmosphere_contributions(argsall...) # all sky
    α_ATM_K, α_K_SFC = surface_atmosphere_contributions(argsclr...) # clear sky

    # This is the normalization done to have energetically consistent albedo value
    # according to Donohoe&Battisti model
    C2_timeavg = α_ATM .- α_ATM_K
    a = C2_timeavg; τ = timemean(Tcorrected); f = timemean(F)
    _x = -2a ./ (τ .* a - f .* τ)
    g = clamp.(1 .- _x ./ Float32(√3), 0, 1)

    C = _cloud_effective_albedo(F, Tcorrected, g)
end

# This actually implements the formula from Lacis
function _cloud_effective_albedo(f, τ, g = 0.9; frequencies = [1.0, 2.0])
    if any(ismissing, τ) # missing values, perform regularization
        @assert hasdim(τ, Time)
        T = sinusoidal_continuation(τ, frequencies)
    else
        T = τ
    end
    p = @. Float32(sqrt(3))*(1 - g)
    R = p .* T ./ (2 .+ p .* T)
    cloud_frac = any(x -> x > 1, f) ? f ./ 100 : f
    eca = cloud_frac .* R
    gstr = g isa Real ? string(g) : "array"
    return ClimArray(Float32.(eca.data), eca.dims, "effective cloud albedo, g = $(gstr)")
end

function cloud_albedo_lacis(f, τ, g)
    p = eltype(g)(sqrt(3))*(1 .- g)
    r = p .* τ ./ (2 .+ p .* τ)
    return f .* r
end
