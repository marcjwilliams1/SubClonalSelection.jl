#Run analysis with neutral data, input parameters:
# mu = 20.0
# cm = 200
# cellularity = 0.7
println()
println("####################")
println("Checking that inference on neutral simulated data with known input parameters returns ground truth...")
srand(1)
out = fitABCmodels("data/neutral.txt",
  "neutral",
  read_depth = 150,
  resultsdirectory = "example",
  nparticles = 200,
  maxiterations = 3 * 10^4,
  Nmax = 10^3,
  maxclones = 1,
  firstpass = false,
  progress = true);

# check if we get the correct model
println("")
println("\tTesting posterior model probability returns neutral as most probable model...")
@test indmax(out.ModelProb[:Probability]) == 1

# extract parameters
mu = out.Posterior[1].Parameters[:mu]
clonalmutations = out.Posterior[1].Parameters[:clonalmutations]
cellularity = out.Posterior[1].Parameters[:cellularity]

println("\tChecking true parameters are within the 80% credible interval range...")
# check if parameters are within 80% credible interval
@test quantile(mu, 0.1) < 20.0 < quantile(mu, 0.9)
@test quantile(clonalmutations, 0.1) < 200.0 < quantile(clonalmutations, 0.9)
@test quantile(cellularity, 0.1) < 0.7 < quantile(cellularity, 0.9)

println("All tests passed for neutral simulated data.")
