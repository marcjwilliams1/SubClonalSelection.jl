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
using StatsBase

import Base.show, ApproxBayes.show

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
  timefuncrand,
  saveresults,
  show

include("types.jl")
include("util.jl")
include("plotting.jl")
include("runabc.jl")

end # module
