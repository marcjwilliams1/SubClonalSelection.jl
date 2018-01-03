# SubClonalSelection

This section will go through a few detailed examples of how a VAF distribution from bulk sequencing data is fitted using Approximate Bayesian Computation. We will go through 3 examples, 2 simulated datasets where the ground truth is known and a final example using data from Nik-Zainal et al. Cell 2012 which is presented in the paper in figure 3.

First I'll go through some of the basics of the package and the options available.

## `fitABCmodels` function
An analysis is performed with the ```fitABCmodels``` function that takes as input either a vector of Floats or a string pointing to a text file containing a vector of floats. These values are the VAF values from a deep sequencing experiment. There are a number of parameters that can be modified in the ```fitABCmodels``` function, these are given below and can also be accessed with ```?fitABCmodels``` in a Julia session.

...
### Function arguments
- `read_depth = 200.0`: Mean read depth of the target data set
- `minreads = 5`: Minimum number of reads to identify a mutation
- `fmin = 0.01`: Minimum range of VAF to perform inference
- `fmax = 0.75`: Maximum range of VAF to perform inference
- `maxiterations = 10^4 `: Maximum number of iterations before ABC terminates
- `maxclones = 2`: Maximum number of clones, can be 1 or 2
- `nparticles = 500`: Number of particles (ie samples) in the ABC output
- `Nmax = 10^4`: Maximum population size used to fit data, increase if suspect that there is a late arising clone
- `resultsdirectory = "output" `: Directory where posterior distributions will be saved to, this directory will be created if it does not already exist.
- `progress = true`: Show progress of ABC sampler with `ProgressMeter` package
- `verbose = true`: Print out summary at each ABC population step.
- `save = false `: Save output or not
- `ϵ1 = 10^6 `: Target ϵ for first ABC step
- `firstpass = false`: If set to true will run a limited first pass of the algorithm to determine a good starting ϵ1 if this unknown.
- `Nmaxinf = 10^10`: Scales selection coefficient value assuming the tumour size is Nmaxinf. Once value >10^9 has limited effect.
- `scalefactor = 2`: Parameter for perturbation kernel for parameter values. Larger values means space will be explored more slowly but fewer particles will be perturbed outside prior range.
- `ρ = 0.0`: Overdispersion parameter for beta-binomial model of sequencing data. ρ = 0.0 means model is binomial sampling
- `adaptpriors = false`: If true priors on μ and clonalmutations are adapted based on the number of mutations in the data set
- `timefunction = timefunc`: Function for KMC algorithm timestep. timefunc returns 1 meaning the timestep is the average of stochastic process. Alternatively timefuncrand can be specified which uses `-log(rand())` to increase the time step, so it is exponentially distributed rather than the mean of the exponential distribution.
- `ploidy = 2`: ploidy of the genome
- `d = 0.0`: Death rate of the thost population in the tumour
- `b = log(2)`: Birth rate of the population. Default is set to `log(2)` so that tumour doubles with each unit increase in t in the absence of cell death.
...

## Example 1 - Neutral synthetic data
For the first example we'll take some synthetic data generated from a neutral simulation of tumour evolution. The input parameters are as follows:
