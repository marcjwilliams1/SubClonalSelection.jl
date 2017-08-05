type Posterior
  MeanHistogram::DataFrame
  Parameters::DataFrame
  Probability::Float64
end

type Results
  ABCsetup
  ABCresults
  VAF::Array{Float64, 1}
  Posterior::Array{Posterior, 1}
  ModelProb::DataFrame
  Sname::String
end
