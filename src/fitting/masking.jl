"""
    normfield(X)
Normalize the field to approximately [0,1] based on the extrema of its
spatial variability.
"""
function normfield(x)
    mi, ma = extrema((timemean(x)))
    return (x .- mi) ./ (ma - mi)
end


"""
	to_non_nan_lats(A)
Given `A` that must have a latitude dimension, return `B` which has the latitudes of
`A` that *do not* contain `NaN` values. Notice that if `A` has a time dimension,
a time-average is done first to return latitudes that would be non-NaN even after
time averaging.
"""
function to_non_nan_lats(A)
	B = hasdim(A, Ti) ? timemean(A) : A
	is = findall(!isnan, B)
	return A[Lat(is)]
end

function dropnan(x)
	is = findall(!isnan, x)
	x[is]
end

oceany_zonalmean(X, OCEAN_MASK) = to_non_nan_lats(timemean(zonalmean(X, OCEAN_MASK)))

"""
    ocean_masked(X, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG) → masked::Vector
Get the field values masked by the Boolean ocean mask,
and limited within `(MINDEG, MAXDEG)`.
"""
function ocean_masked(X, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG)
    latitude_bounding = Coord(Lat(Between(MINDEG, MAXDEG)))
    Xbounded = X[latitude_bounding]
    bounded_mask_selection = OCEAN_MASK[latitude_bounding]
    return Xbounded.data[bounded_mask_selection.data]
end


"""
    maskedtimezonalmean(X, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG) → masked::Vector
The "zonal and time mean of `X` over ocean",
and limited within `(MINDEG, MAXDEG)`.
"""
function maskedtimezonalmean(Φ, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG)
    to_non_nan_lats(
        timemean(
            zonalmean(
                Φ[Coord(Lat((MINDEG)..(MAXDEG)))],
                OCEAN_MASK[Coord(Lat((MINDEG)..(MAXDEG)))]
            )
        )
    )
end


"""
    maskedspacemean(X, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG) → masked::Vector
Return the spatial map (i.e., timemean), masked by the given mask,
and limited within `(MINDEG, MAXDEG)`.
"""
function maskedspacemean(Φ, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG)
    x = timemean(
        Φ[Coord(Lat(MINDEG..MAXDEG))],
        OCEAN_MASK[Coord(Lat(MINDEG..MAXDEG))]
    )
    return dropnan(x)
end


"""
    maskedtimezonalmean(X, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG) → masked::Vector
Given a timemean of a field do a masked zonal mean,
and limited within `(MINDEG, MAXDEG)`.
"""
function maskedzonalmean(Φ, OCEAN_MASK, MAXDEG, MINDEG = -MAXDEG)
    to_non_nan_lats(
        zonalmean(
            Φ[Coord(Lat(MINDEG..MAXDEG))],
            timemean(OCEAN_MASK[Coord(Lat(MINDEG..MAXDEG))])
        )
    )
end
