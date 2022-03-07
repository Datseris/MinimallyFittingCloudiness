using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
include(srcdir("fitting", "general.jl"))

predictors = (:Ω_mean, :ECTEI)
predicted = :C
Φname = "\$C\$"
expression = "p[1]*50*(tanh(p[2]*x1 + p[3]*x2) + 1)"
params = [0.40704233122044947, 6.870173231618394, 0.0803937826587933]
MAXDEG = 70 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
limits = (0, 40)
seasonal_limits = (-4, 14)
seasonal_offset = 11

function generate_contributions(p, Ps)
    c1 = p[2]*Ps[1]
    c2 = p[3]*Ps[2]
    c3 = p[2]*Ps[1] .+ p[3]*Ps[2]
    c1 = @. p[1]*50*(tanh(p[2]*Ps[1]) + 1)
    c2 = @. p[1]*50*(tanh(p[3]*Ps[2]) + 1)
    return c1, c2, c3
end
# contribution_titles = ("\$(\\tanh(p_2 \\Omega)+1)/2\$", "\$p_1(\\tanh(p_2 \\Omega)+1)/2\$", "\$p_2 I\$")
contribution_titles = ("\$50 p_1(\\tanh(p_2\\Omega) + 1)\$", "\$50 p_1(\\tanh(p_3 I) + 1)\$", "\$p_2 \\Omega + p_3 I\$")
contrib_cmaps = (:inferno, :inferno, :BrBG)

# The rest of this will be made on script
include(papersdir("best_fit_plot.jl"))