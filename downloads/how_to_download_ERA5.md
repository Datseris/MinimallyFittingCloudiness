To download data from Copernicus we use PyCall.jl and its Python download API.
The API is Described here: https://cds.climate.copernicus.eu/api-how-to
Notice that you need to set up an API key as instructed in the above webpage,
(i.e., create a document in `$HOME/.cdsapirc` that contains your key).

While when selecting data online ("Datasets" tab), there is a button at
the end that says "Show API request". This is based on the Python package `cdsapi`.
We install it in Julia using PyCall.jl:
```julia
using PyCall
PyCall.Conda.add("cdsapi"; channel = "conda-forge")
```
That's it. The provided Julia scripts use transform Julia dictionaries to 
Python dictionaries and call the API. They are super convenient. 

Just be sure to set properly the `datapath` variable, which is the path to
the folder to save the ERA5 data. Use ClimateBase.jl to open the data!
