using ClimateBase, Dates
# This needs to numpy, xarray and dask installed from Conda
using PyCall
xr = pyimport("xarray")
np = pyimport("numpy")

"""
    climarray_from_xarray(xa, fieldname, name = Symbol(fieldname))
Load underlying field with given `fieldname` from the given xarray instance `xa`,
optionally providing a name for it. This `xa` can be loaded with commands like
```julia
using PyCall, ClimateBase
xr = pyimport("xarray")
ERA5_files = "some_file_name.nc"
xa = xr.open_mfdataset(ERA5_files)
X = climarray_from_xarray(xa, "w", "optional name")
```
"""
function climarray_from_xarray(xa, fieldname, name = Symbol(fieldname))
    w = getproperty(xa, Symbol(fieldname))
    raw_data = Array(w.values)
    dnames = collect(w.dims) # dimensions in string name
    dim_values, dim_attrs = extract_dimension_values_xarray(xa, dnames)
    @assert collect(size(raw_data)) == length.(dim_values)
    actual_dims = create_dims_xarray(dnames, dim_values, dim_attrs)
    ca = ClimArray(raw_data, actual_dims, name; attrib = w.attrs)
end

function extract_dimension_values_xarray(xa, dnames = collect(xa.dims))
    dim_values = []
    dim_attrs = Vector{Any}(fill(nothing, length(dnames)))
    for (i, d) in enumerate(dnames)
        dim_attrs[i] = getproperty(xa, d).attrs
        x = getproperty(xa, d).values
        if d ≠ "time"
            push!(dim_values, x)
        else
            dates = [np.datetime_as_string(y)[1:19] for y in x]
            dates = DateTime.(dates)
            push!(dim_values, dates)
        end
    end
    return dim_values, dim_attrs
end

function create_dims_xarray(dnames, dim_values, dim_attrs)
    true_dims = ClimateBase.to_proper_dimensions(dnames)
    optimal_values = ClimateBase.vector2range.(dim_values)
    out = []
    for i in 1:length(true_dims)
        # push!(out, true_dims[i](optimal_values[i]; metadata = dim_attrs[i]))
        push!(out, true_dims[i](optimal_values[i]))
    end
    return (out...,)
end
