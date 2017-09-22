module SubClonalSelection

# package code goes here
using CancerSeqSim
using DataFrames
using ApproxBayes
using Distributions
using Distances
using GLM
using Gadfly
using Colors

import Base.show

export
  # types
  Results,


  #functions
  fitABCmodels,
  plothistogram,
  plotmodelposterior,
  plotparameterposterior,
  saveallplots,
  timefunc,
  timefuncrand

include("types.jl")
include("runabc.jl")
include("util.jl")
include("plotting.jl")

end # module
