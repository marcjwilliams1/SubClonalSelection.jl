#Run analysis with non-neutral data with 1 clone, input parameters:
# mu = 20.0
# cm = 200
# cellularity = 0.7
println()
println("####################")
println("Checking inference on non-neutral simulated data with known input parameters returns ground truth...")
srand(1)
out = fitABCmodels("data/oneclone.txt",
  "oneclone",
  read_depth = 150,
  resultsdirectory = "example",
  nparticles = 200,
  maxiterations = 3 * 10^4,
  Nmax = 10^3,
  maxclones = 2,
  verbose = true,
  progress = true);

# check if we get the correct model
println("\tTesting posterior model probability returns neutral as most probable model...")
@test out.ModelProb[:Probability][2] > 0.5

# extract parameters
mu = out.Posterior[2].Parameters[:mu]
clonalmutations = out.Posterior[2].Parameters[:clonalmutations]
cellularity = out.Posterior[2].Parameters[:cellularity]

println("\tChecking true parameters are within the 80% credible interval range...")
# check if parameters are within 80% credible interval
@test quantile(mu, 0.1) < 20.0 < quantile(mu, 0.9)
@test quantile(clonalmutations, 0.1) < 200.0 < quantile(clonalmutations, 0.9)
@test quantile(cellularity, 0.1) < 0.8 < quantile(cellularity, 0.9)
