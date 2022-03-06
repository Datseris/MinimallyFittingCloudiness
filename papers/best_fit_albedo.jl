using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
include(srcdir("fitting", "general.jl"))

predictors = (:Ω_nf, :ECTEI)
predicted = :C
Φname = "\$C\$"
expression = "p[1]*x1 + p[2]*x2*(1 - x1)"
params = [49.0908918702647, 2.1032725951879363]
MAXDEG = 70 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
limits = (0, 40)
seasonal_limits = (-4, 14)
seasonal_offset = 11

function generate_contributions(p, Ps)
    c1 = p[1]*Ps[1]
    c2 = p[2]*Ps[2]
    # c3 =  - p[2] .* Ps[2] .* Ps[1]
    c3 =  p[2] .* Ps[2] .* ( 1 .- Ps[1] )
    return c1, c2, c3
end
contribution_titles = ("\$p_1 U\$", "\$p_2 I\$", "\$p_2 I (1-U)\$")

# The rest of this will be made on script
include(papersdir("best_fit_plot.jl"))