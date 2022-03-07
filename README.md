# MinimallyFittingCloudiness

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named accompanying the paper _Minimalistic fits to planetary cloudiness_, by Datseris et al.


**TOC**
1. [Reproducing the code](#reproducing-the-code)
2. [Using the code base](#using-the-code-base)

## Reproducing the code
To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts. Notice that plotting of maps is done using the `cartopy` package of Python. To this end, you would need to run the commands:
```julia
Pkg.add("Conda")
import Conda
Conda.add("cartopy")
```

## Using the code base
This code base provides the means to fit arbitrary functions of arbitrary predictors versus a field-to-be-predicted, using the method described in our paper (Section 2.3). The source code and a premade script make this possible. The script `scripts/field_definitions.jl` defines all spatiotemporal fields in a form of a dictionary. The script `scripts/general_model_fit.jl` allows the user to define and run arbitrary models while only changing three lines of code.

**It is crucial that all data are in an equal area grid.** Using CDO, this can be done with the command `cdo remapbil,gea250 -setgridtype,lonlat IN.nc OUT.nc`. It is probably much simpler if you just send me an email asking for the data already mapped into an equal area grid.

The code that actually makes the figures of the paper in is `papers`. The code does not use the same naming convention as the paper.
