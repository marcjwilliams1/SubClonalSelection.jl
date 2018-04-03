# SubClonalSelection

This section will go through a few detailed examples of how a VAF distribution from bulk sequencing data is fitted using Approximate Bayesian Computation. We will go through 3 examples, 2 simulated datasets where the ground truth is known and a final example using data from Zhang et al. Science 2015 which is presented in the paper in figure 3.

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
- `adaptpriors = false`: If true priors on μ and clonalmutations are adapted based on the number of mutations in the data set. This is an experimental feature that needs further validation, although initial tests suggest it performs well and does not skew inferences. To run the inference with default priors as described in the paper keep this to false.
- `timefunction = timefunc`: Function for KMC algorithm timestep. timefunc returns 1 meaning the timestep is the average of stochastic process. Alternatively timefuncrand can be specified which uses `-log(rand())` to increase the time step, so it is exponentially distributed rather than the mean of the exponential distribution.
- `ploidy = 2`: ploidy of the genome
- `d = 0.0`: Death rate of the host population in the tumour, we would advise keeping this to 0.0 as it is only μ/β that can be inferred. That is differences in the death rate are unidentifiable in this framework.
- `b = log(2)`: Birth rate of the population. Default is set to `log(2)` so that tumour doubles with each unit increase in t in the absence of cell death.

## Example 1 - Neutral synthetic data
For the first example we'll take some synthetic data ("neutral.txt") generated from a neutral simulation of tumour evolution. The input parameters for the simulations were as follows:
- Mutation rate: 20.0
- Number of clonal mutation: 300
- Number of subclones: 0 (ie neutral)
- Cellularity: 0.7
- Tumour population size: 10^6
- Read depth: 300X

First we'll load the packages that we need.
```julia
using SubClonalSelection
using Gadfly
using DataFrames
```
We'll now use ```fitABCmodels``` from ```SubClonalSelection``` to attempt to recover these parameters as well as the number of subclones.

```julia
srand(123)
out = fitABCmodels("example/neutral.txt",
    "neutral",
    read_depth = 300,
    resultsdirectory = "example/",
    nparticles = 100,
    maxiterations = 4*10^4,
    maxclones = 1,
    save = true,
    adaptpriors = true,
    verbose = true,
    Nmaxinf = 10^6,
    fmin = 0.01)
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
- Mutations in subclone: 251
- Time emerges (tumour doublings): 9.0
- Read depth: 300X

As before we'll use ```fitABCmodels``` to attempt to recover these parameters as well as the correct number of subclones (1).

```julia
srand(123)
@time out1 = fitABCmodels("example/oneclone.txt",
  "oneclone",
  read_depth = 300,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 4*10^4,
  maxclones = 2,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^6,
  fmin = 0.01);
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

## Example 3 - Lung cancer sample
This is an example from data we presented in the paper in figure 3C, for this sample we found evidence for one subclone. We measured the overdispersion parameter of this data to be 0.005 which we can input into the inference algorithm. This data has also been corrected for the purity of the sample so we constrain this in our inference (it is still adviseable to give this a bit of freedom, hence we set the lower limit to 0.95). We also notice the mode of the lower peak ~ 0.04 so we set minvaf to this number.

```julia
srand(123)
@time out = fitABCmodels("example/4990-12/data/4990-12.txt",
  "4990-12",
  read_depth = 150,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 4*10^4,
  maxclones = 2,
  ρ = 0.005, #measured overdispersion using clonal peak
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^10,
  minvaf = 0.04, #minimum vaf to resolve mutations, can be seen from lower peak
  mincellularity = 0.95);
```

First we can see that our model can accurately fit the data.
```julia
plothistogram(out, 1) #1 specified to only plot data from simulations of model 0
```
![plot](/example/4990-12/plots/4990-12-histogram-1clone.png)

And secondly that the 1 clone model is favoured as we found in the paper.
```julia
plotmodelposterior(out)
```
![plot](/example/4990-12/plots/4990-12-modelposterior.png)

Finally we can plot the inferred parameters from the model.
```julia
plotparameterposterior(out, 1)
```
![plot](/example/4990-12/plots/4990-12-posterior-1clone.png)

### Notes
We note that in these examples above we have used a limited number of particles and a limited number of iterations to what would normally be recommended and what was used in the paper, we used these examples as it should be feasible to run them on a laptop in <30 minutes or so and as should be apparent the results are reasonably good. Nonetheless in general we would recommend running the algorithm with 500 particles/samples and for a minimum of 10^6 iterations which is computationally expensive and hence is best suited to a HPC of some sort. Note that the inferences improves with increasing the number of iterations as the error between the target data set and the simulated datasets decreases. This is particularly relevant when considering up to 2 subclones (here the search space is large) and a large number of iterations is required to correctly identify samples with 2 subclones and their corresponding parameters.

Also note that the reported mutation rate is the effective mutation rate per tumour doubling. If you want to convert this to a quantity in terms of the per bp per tumour doubling you'll need to divide this number by the size of the target that was sequenced (eg ~ 30*10^6 for WXS).

There are plans to parallelise the algorithm in future which should help the speed.
