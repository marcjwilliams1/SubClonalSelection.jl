# SubClonalSelection

This section will go through a few detailed examples of how a VAF distribution from bulk sequencing data is fitted using Approximate Bayesian Computation. We will go through 3 examples, 2 simulated datasets where the ground truth is known and a final example using data from Nik-Zainal et al. Cell 2012 which is presented in the paper in figure 3.

First I'll go through some of the basics of the package and the options available.

## `fitABCmodels` function
An analysis is performed with the ```fitABCmodels``` function that takes as input either a vector of Floats or a string pointing to a text file containing a vector of floats. These values are the VAF values from a deep sequencing experiment. There are a number of parameters that can be modified in the ```fitABCmodels``` function, these are given below and can also be accessed with ```?fitABCmodels``` in a Julia session.


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

## Example 1 - Neutral synthetic data
For the first example we'll take some synthetic data ("neutral.txt") generated from a neutral simulation of tumour evolution. The input parameters for the simulations were as follows:
- Mutation rate: 20.0
- Number of clonal mutation: 300
- Number of subclones: 0 (ie neutral)
- Cellularity: 0.7
- Tumour population size: 10^6

First we'll load the packages that we need.
```julia
using SubClonalSelection
using Gadfly
using DataFrames
```
We'll now use ```fitABCmodels``` from ```SubClonalSelection``` to attempt to recover these parameters as well as the number of subclones.

```julia
out = fitABCmodels("example/neutral.txt", # text file with data
  "neutral", # sample name
  read_depth = 150,
  resultsdirectory = "", #use this directory
  nparticles = 100,
  maxiterations = 2 * 10^5,
  Nmax = 10^3,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true,
  Nmaxinf = 10^6);
```

This may take ~30-60 minutes on a desktop computer. With this output we can then plot the distributions to see if we get the right answers. First we'll plot the model posterior probabilities, we would hope to see model 0 (0 subclones) with the highest probability which is exactly what we see.
```julia
plotmodelposterior(out)
```
![plot](/example/neutral/plots/neutral-modelposterior.png)

We can visually inspect how well this model fits the data by overlaying a summary (mean and 95% credible intervals) of the VAF data from simulations that were accepted on top of the target data set.
Plot the histogram for model 2.
```julia
plothistogram(out, 0) #0 specified to only plot data from simulations of model 0
```
![plot](/example/neutral/plots/neutral-histogram-0clone.png)
Finally we can plot the posterior distributions of the parameters and check whether we have correctly identified the true parameters. As would be hoped in this case the mode of the posteriors do indeed closely match the true input parameters.
```julia
plotparameterposterior(out, 0)
```
![plot](/example/neutral/plots/neutral-posterior-0clone.png)


## Example 2 - Synthetic data with 1 subclone
For this second example  we'll take some synthetic data ("oneclone.txt") generated from a simulation of tumour evolution that contains a single subclone. The input parameters for the simulations were as follows:
- Mutation rate: 20.0
- Number of clonal mutation: 300
- Number of subclones: 1
- Cellularity: 0.7
- Tumour population size: 10^6
- Subclone frequency: 0.58
- Fitness advantage: 1.03
- Time emerges (tumour doublings): 9.0

As before we'll use ```fitABCmodels``` to attempt to recover these parameters as well as the correct number of subclones (1).

```julia
out = fitABCmodels("example/oneclone.txt", # text file with data
  "oneclone", # sample name
  read_depth = 150,
  resultsdirectory = "", #use this directory
  nparticles = 100,
  maxiterations = 2 * 10^5,
  Nmax = 10^3,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true,
  Nmaxinf = 10^6);
```

As in the neutral example we can confirm we recover the input number of subclones and parameters by plotting the posterior distribution.
```julia
plotmodelposterior(out)
```
![plot](/example/oneclone/plots/oneclone-modelposterior.png)

```julia
plothistogram(out, 1) #1 specified to only plot data from simulations of model 0
```
![plot](/example/oneclone/plots/oneclone-histogram-1clone.png)

```julia
plotparameterposterior(out, 1)
```
![plot](/example/oneclone/plots/oneclone-posterior-1clone.png)


## Example 3
For the final example we'll take some real data from Nik-Zainal et al Cell 2012 which is presented in figure 3 Williams et al 2018. We found that this data was most consistent with 2 subclones, which can be easily seen from the data.

```julia
out = fitABCmodels("example/nikzainal.txt",
  "nikzainal",
  read_depth = 180,
  resultsdirectory = "",
  nparticles = 100,
  maxiterations = 2 * 10^5,
  minreads = 7, # we observed that we could detect things above ~4% VAF
  Nmax = 10^3,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true);
```

We'll now look at the results of the inference as before.
```julia
plotmodelposterior(out)
```
![plot](/example/nikzainal/plots/nikzainal-modelposterior.png)

```julia
plothistogram(out, 2) #2 specified to only plot data from simulations of model 0
```
![plot](/example/nikzainal/plots/nikzainal-histogram-2clone.png)

```julia
plotparameterposterior(out, 2)
```
![plot](/example/nikzainal/plots/nikzainal-posterior-2clone.png)
