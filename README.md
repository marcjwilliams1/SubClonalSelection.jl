# SubClonalSelection

A Julia package for inferring the strength of selection from cancer sequencing data. Package simultaneously estimates the number of subclones in the population and their relative fitnesses. This is done by generating synthetic data by simulating different population dynamics (see [CancerSeqSim.jl](https://github.com/marcjwilliams1/ApproxBayes.jl)) and fitting this to the model using Approximate Bayesian Computation (see [ApproxBayes.jl](https://github.com/marcjwilliams1/CancerSeqSim.jl)).

## Getting Started
Package has been tested extensively with [Julia](https://julialang.org/) v0.5.1 but should work with later versions. If there any problems please report an issue.

To download the package, once you're in a Julia session type the following command:
```
Pkg.clone("https://github.com/marcjwilliams1/SubClonalSelection.jl")
```

You will also need the [CancerSeqSim.jl](https://github.com/marcjwilliams1/ApproxBayes.jl) and  [ApproxBayes.jl](https://github.com/marcjwilliams1/CancerSeqSim.jl) packages, which can be downloaded using the following (all other dependencies should be downloaded automatically).
```
Pkg.clone("https://github.com/marcjwilliams1/ApproxBayes.jl")
Pkg.clone("https://github.com/marcjwilliams1/CancerSeqSim.jl")
```

## Input data
Running an analysis requires variant allele frequencies (VAFs) as measured in deep sequencing of cancer samples. Generating the synthetic data assumes that the cancer is diploid, therefore any mutations falling in non-diploid regions should be removed. This does unfortunately mean that in highly aneuploid tumours, there will not be enough mutations to perform an analysis. We would recommend a minimum of 100 mutations.

## Running an analysis
The main function to perform an analysis is ```fitABCmodels```. This takes as its first argument either a vector of Floats or a string pointing to a text file with a vector of floats, which will be read in automatically. The second argument is the name of the sample you wish to analyse which will be used to write the data and plots to a file. There are then a number of keyword arguments set to reasonable defaults. More details of these arguments and their defaults can be found by typing ```?fitABCmodels``` in a Julia session.

There is some example data generated from the simulation found in the examples directory. For example the following command will run the inference on the ```oneclone.txt``` data set with 200 posterior samples and 5*10^4 iterations given sequencing depth of sample is 150X:
```
out = fitABCmodels("example/oneclone.txt",
  "oneclone",
  read_depth = 150,
  resultsdirectory = "example",
  nparticles = 200,
  maxiterations = 5*10^4);
```
The above command should run in about 15 minutes on a reasonably specced computer. For robust inferences we would recommend using 500 particles and setting ```maxiterations``` to 10^6. This starts to get computationally expensive so running on a cluster would be recommended in most cases.

Also included are a number of functions to summarize the output and plot the posterior. ```show(out)``` will print a summary of the posterior model and parameter probabilities. We can also plot the posterior distributions.

Plot the posterior model probabilities.
```
plotmodelposterior(out)
```
![plot](/example/oneclone/plots/modelposterior.png)

Plot the histogram for model 2.
```
plothistogram(out, 1)
```
![plot](/example/oneclone/plots/histogram-1clone.png)

Plot the posterior parameter distribution for model 2.
```
plothistogram(out, 1)
```
![plot](/example/oneclone/plots/posterior-1clone.png)

Note the ground truth of the parameters in this case are:
```
  mu = 20.0
  clonalmutations = 200
  s = 0.6
  t = 8.3
  cellularity = 0.8
  freq = 0.59
  scmuts = 180
```

Finally we can also save all plots and text files with posterior distributions to a directory, unless specified the default will be to create a file a directory called ```output``` in the current working directory.

```
saveresults(out; resultsdirectory = "example")
saveallplots(out, resultsdirectory = "example")
```

### Additional info
One issue with model selection in the ABC SMC framework is that sometimes models can die out. This is particularly an issue if the prior probability of a particular parameter set is low, this will be the case for models with many parameters due to the large parameter space. For this application this can be an issue for the 2 clone model, which has the largest parameter space. To mitigate the problem of this model dying out before the algorithm has converged (where the most likely model could well be the 2 clone model) there is the optional argument `firstpass`, which if set to `true` will run a first pass on the data to choose a reasonable starting ϵ, where parameters already fit the data to a reasonable level before the SMC algorithm commences. This mitigates the issue somewhat by in effect reducing the search space. Despite this, sometimes the 2 clone model will still die out, so for this reason we would recommend running the algorithm multiple times (this is recommended in any case and is good practice in Bayesian computational approaches).

Note that if there is strong evidence for a neutral model, then having the 2 clone model with 0 probability is not a problem and in fact is quite likely. The additional model complexity provided by the selection model is unnecessary to explain the data.
