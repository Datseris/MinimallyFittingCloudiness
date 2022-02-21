# MinimallyFittingCloudiness

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> MinimallyFittingCloudiness

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

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Using the code base
This code base provides the source code and a premade script (or Jupyter notebook) that allows one to fit any arbitrary spatiotemporal field with any arbitrary function of any other arbitrary spatiotemporal fields, as we have done in our paper.

The current status of the scripts of course are for fitting cloudiness, but one can change the predictors or fields to be predicted by adjusting the `scripts/fields_definitions.jl` file, and then just updating the fields you want to use in the `scripts/general_model_fit.jl` file. 

**ABSOLUTELY CRUCIAL FOR FIELDS TO BE IN EQUAL AREA GRID**