# SubClonalSelection

[![Build Status](https://travis-ci.org/marcjwilliams1/SubClonalSelection.jl.svg?branch=master)](https://travis-ci.org/marcjwilliams1/SubClonalSelection.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/marcjwilliams1/SubClonalSelection.jl?branch=master&svg=true)](https://ci.appveyor.com/project/marcjwilliams1/SubClonalSelection-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/marcjwilliams1/SubClonalSelection.jl/badge.svg?branch=master)](https://coveralls.io/github/marcjwilliams1/SubClonalSelection.jl?branch=master)
[![codecov.io](http://codecov.io/github/marcjwilliams1/SubClonalSelection.jl/coverage.svg?branch=master)](http://codecov.io/github/marcjwilliams1/SubClonalSelection.jl?branch=master)

A Julia package for inferring the strength of selection from cancer sequencing data. Package simultaneously estimates the number of subclones in the population and their relative fitnesses. This is done by generating synthetic data by simulating different population dynamics (see [CancerSeqSim.jl](https://github.com/marcjwilliams1/CancerSeqSim.jl)) and fitting this to the model using Approximate Bayesian Computation (see [ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl)).

## Getting Started
Package has been tested extensively with [Julia](https://julialang.org/) v0.5.1 but should work with later versions. If there any problems please report an issue.

To download the package, once you're in a Julia session type the following command:
```
Pkg.clone("https://github.com/marcjwilliams1/SubClonalSelection.jl")
```

You will also need the [ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl) and  [CancerSeqSim.jl](https://github.com/marcjwilliams1/CancerSeqSim.jl) packages, which can be downloaded using the following (all other dependencies should be downloaded automatically).
```
Pkg.clone("https://github.com/marcjwilliams1/ApproxBayes.jl")
Pkg.clone("https://github.com/marcjwilliams1/CancerSeqSim.jl")
```

Then once you are in a julia session the package can be loaded with
```
using SubClonalSelection
```
Running `Pkg.test("SubClonalSelection") will run a test suite to confirm the algorithm works on some test data and recovers the ground truth from some know synthetic data.

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
