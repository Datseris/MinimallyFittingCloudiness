using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(srcdir("plotting", "style.jl"))

error2_name = "mean_pearson"  # upper-left corner
error1_name = "full_error" # lower-right corner
error1_name = "timezonal_error" # lower-right corner
error2_name = "seasonal_error"  # upper-left corner

fig = figure()
axwidth = 20
allow_intercept = false

names_to_symbols = Dict(
    "Ω_mean" => "\\Omega",
    "Ω_std" => "S",
    "Ω_nf" => "U",
    "WS10" => "W",
    "Tsfc" => "T",
    "EIS" => "E",
    "ECTEI" => "I",
    "V" => "V",
    "q700" => "Q",
)

titlenames = Dict("C"=> "cloud albedo, \$C\$", "L"=>"cloud longwave rad. eff., \$L\$")
pcm = nothing
textsize = 22


for (axi, predicted) in enumerate((:C, :L))

d = wload(datadir("modelfits", "linear_2_predictors", "$(predicted)_intercept=$(allow_intercept).jld2"))
@unpack predictors, = d
pnames = ["\$"*names_to_symbols[string(p)]*"\$" for p in predictors]

function highlight_cell(i; ax=gca(), kwargs...)
    x,y = i
    rect = plt.Rectangle((x-.5, y-.5), 1,1; fill=false, color = "black", lw=3, kwargs...)
    ax.add_patch(rect)
end

function mink(a, k)
    linind = CartesianIndices(a)
    b = partialsortperm(vec(a), 1:k)
    return [linind[x] for x in b]
end

error1 = d[error1_name]
error2 = d[error2_name]

toplot = copy(error1)
k = 3
min1 = mink(error1, k)
min2 = mink(error2, k)
L = size(error1, 1)

for i in 1:L
    for j in (i+1):L
        toplot[j,i] = error2[i,j]
    end
    toplot[i,i] = NaN
end

ax = subplot2grid((1, 2axwidth + 1), (0, (axi-1)*axwidth); colspan = axwidth)
x = (1:length(pnames)+1) .- 0.5
pcm = ax.pcolormesh(x, x, toplot; vmin = 0, vmax = 1, cmap=plt.cm.get_cmap("viridis", 10))
ax.set_title(titlenames[string(predicted)])
xticks(x[1:end-1] .+ 0.5, pnames; size = textsize)
if axi == 1
    yticks(x[1:end-1] .+ 0.5, pnames; size = textsize)
else
    yticks(x[1:end-1] .+ 0.5, [])
end

for i in 1:k
    # Becuse puthon arrays are plotted as transposed, i need to reverse the
    # opposite indices... # TODO: Switch to makie.
    highlight_cell(reverse(min1[i].I); color = :red)
    highlight_cell(min2[i].I; color = :red)
end
j = argmin(error1 .* error2)
highlight_cell(reverse(j.I); color = :black, linewidth = 4, linestyle = :dashed)
highlight_cell(j.I; color = :black, linewidth = 4, linestyle = :dashed)

# if axi == 1
#     text(0.95, 0.05, error1_name; transform = ax.transAxes, size = textsize, ha="right")
#     text(0.05, 0.95, error2_name; transform = ax.transAxes, size = textsize)
# end

end


axcb = subplot2grid((1, 2axwidth + 1), (0, 2*axwidth))
colorbar(pcm; label = "error, \$\\epsilon\$", cax = axcb, extend = "max")
fig.subplots_adjust(wspace = 0.8, left = 0.05, right = 0.9)
wsave(papersdir("plots", "two_predictors_linear"), fig)