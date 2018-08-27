mutable struct Posterior
  MeanHistogram::DataFrame
  Parameters::DataFrame
  Probability::Float64
end

mutable struct Results
  ABCsetup
  ABCresults
  VAF::Array{Float64, 1}
  Posterior::Array{Posterior, 1}
  ModelProb::DataFrame
  SampleName::String
end
