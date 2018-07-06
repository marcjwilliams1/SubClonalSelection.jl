#Run analysis with non-neutral data with 1 clone, input parameters:
# mu = 20.0
# cm = 200
# cellularity = 0.8
println()
println("####################")
println("Checking inference on non-neutral simulated data with known input parameters returns ground truth...")
srand(1)
@time out = fitABCmodels("data/oneclone.txt",
  "oneclone",
  read_depth = 150,
  resultsdirectory = "output",
  nparticles = 100,
  maxiterations =  10^4,
  Nmax = 10^3,
  maxclones = 2,
  verbose = true,
  firstpass = false,
  save = true,
  progress = true);

# check if we get the correct model
println("\tTesting posterior model probability returns one clone as most probable model...")
@test indmax(out.ModelProb[:Probability]) == 2

println("\tTesting 2 clone model did not die out...")
@test out.ModelProb[:Probability][3] > 0.0

# extract parameters
mu = out.Posterior[2].Parameters[:mu]
clonalmutations = out.Posterior[2].Parameters[:clonalmutations]
cellularity = out.Posterior[2].Parameters[:cellularity]

println("\tChecking true parameters are within the 95% credible interval range...")
# check if parameters are within 80% credible interval
@test quantile(mu, 0.025) < 20.0 < quantile(mu, 0.975)
@test quantile(clonalmutations, 0.025) < 200.0 < quantile(clonalmutations, 0.975)
@test quantile(cellularity, 0.025) < 0.8 < quantile(cellularity, 0.975)

println("\tChecking plotting functions work and are saved...")
saveallplots(out)
@test isfile("output/oneclone/plots/oneclone-histogram-0clone.pdf")
@test isfile("output/oneclone/posterior/oneclone-histogram-clone0.csv")
rm("output", recursive = true)


println("All tests passed for single clone.")
