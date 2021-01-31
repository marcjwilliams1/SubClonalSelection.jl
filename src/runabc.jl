function gettargetDF(VAF; fmin = 0.05, fmax = 0.75)
  targetdataCDF = getCDF(VAF, 0.001, fmin = fmin, fmax = fmax)
  return targetdataCDF, VAF
end

function getmodel(abcres, model)
    indeces = map(x -> x.model, abcres.particles) .== model
    abcresnew = deepcopy(abcres)
    abcresnew.particles = abcresnew.particles[indeces]
    abcresnew.parameters = abcresnew.parameters[model]
    abcresnew.weights = abcresnew.weights[model]
    return abcresnew
end

function getCDF(VAF::Array, step_size::Float64; fmin = 0.05, fmax = 0.7)
  #fast way to calculate CDF
  out = cumsum(fit(Histogram, VAF, fmax:-step_size:fmin,closed=:left).weights[1:end - 1])
  out = out .- out[1]
  return out
end

function tumourABCneutral(parameters, constants, targetdata)
    #function to simulate neutral tumour growth, outputs the euclidean distance between targetdata and simulation and
    #and histogram values which can be used to plot the output

    cst = constants
    simdata = simulate(nclones = 0,
                ploidy = cst[1],
                read_depth = cst[2],
                μ = parameters[1],
                clonalmutations = round(parameters[2]),
                d = cst[3],
                b = cst[4],
                ρ = cst[5],
                Nmax = 500,
                timefunction = cst[7],
                detectionlimit = cst[8],
                s = Float64[],
                tevent = Float64[],
                cellularity = parameters[3])

    simdatacdf = getCDF(simdata.sampleddata.VAF, 0.001,
                    fmin = cst[9],
                    fmax = cst[10])

    return euclidean(simdatacdf, targetdata), [simdata.sampleddata.DF], true
end


function tumourABCselection(parameters, constants, targetdata)

    # we only consider the tumour to have a clone if it is observeable in the data, above frequency 5%.

    cst = constants
    simdata = simulate(nclones = 1,
                ploidy = cst[1],
                read_depth = cst[2],
                μ = parameters[1],
                clonalmutations = round(parameters[2]),
                d = cst[3],
                b = cst[4],
                ρ = cst[5],
                Nmax = cst[6],
                s = [parameters[3]],
                tevent = [parameters[4]],
                timefunction = cst[7],
                detectionlimit = cst[8],
                cellularity = parameters[5])

    simdatacdf = getCDF(simdata.sampleddata.VAF, 0.001,
                    fmin = cst[9],
                    fmax = cst[10])

    if length(simdata.output.subclonemutations) == 0
        out = [simdata.sampleddata.DF, 0.0, 0.0, 0.0]
    else
        out = [simdata.sampleddata.DF, simdata.output.subclonemutations[1], simdata.output.Ndivisions[1], simdata.output.clonefreq[1]]
    end

    #return true if clone is between 0.05 and 0.95
    c = false
    if ((sum(simdata.output.clonefreq.<0.95).==1) & (sum(simdata.output.clonefreq.>0.05).==1))[1] == true
        c = true
    end

    euclidean(simdatacdf, targetdata), out, c
end

function tumourABCselection2(parameters, constants, targetdata)

    cst = constants
    simdata = simulate(nclones = 2,
                ploidy = cst[1],
                read_depth = cst[2],
                μ = parameters[1],
                clonalmutations = round(parameters[2]),
                d = cst[3],
                b = cst[4],
                ρ = cst[5],
                Nmax = cst[6],
                s = [parameters[3], parameters[5]],
                tevent = [parameters[4], parameters[6]],
                timefunction = cst[7],
                detectionlimit = cst[8],
                cellularity = parameters[7])

    simdatacdf = getCDF(simdata.sampleddata.VAF, 0.001,
                    fmin = cst[9],
                    fmax = cst[10])

    if length(simdata.output.subclonemutations) < 2
        out = [simdata.sampleddata.DF, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    else
        out = [simdata.sampleddata.DF, simdata.output.subclonemutations[1],
        simdata.output.subclonemutations[2],
        simdata.output.Ndivisions[1],
        simdata.output.Ndivisions[2],
        simdata.output.clonefreq[1],
        simdata.output.clonefreq[2]]
    end

    #return true if clones are between 0.05 and 0.95
    c = false
    if ((sum(simdata.output.clonefreq.<0.95).==2) & (sum(simdata.output.clonefreq.>0.05).==2))[1] == true
        c = true
    end

    euclidean(simdatacdf, targetdata), out, c
end

timefunc() = 1
timefuncrand() = -log(rand())

function getsetup(maxclones; nparticles = 100, maxiterations = 10^4, convergence = 0.1, ϵT = 1.0, read_depth = 100.0, Nmax = 10^3, detectionlimit = 0.05, modelkern = 0.4, ϵ1 = 10^6, mincellularity = 0.1, maxcellularity = 1.1, ρ = 0.0, maxclonalmutations = 10000.0, maxmu = 620.0, timefunction = timefunc, ploidy = 2, d = 0.0, b = log(2), fmin = 0.05, fmax = 0.05, Nmaxinf = 10^10)

  #function to define priors, constants and create ABCSMC model type

  read_depth = read_depth
  Nmax = Nmax
  ρ = ρ

  cst = [ploidy, read_depth, d, b, ρ, Nmax, timefunction, detectionlimit, fmin, fmax];

  priormu = [0.01, maxmu]
  priorcm = [0.0, Float64(maxclonalmutations)]
  priorcellularity = [mincellularity, maxcellularity]

  #need to create Prior type which has a distribution type array with a corresponding distribution specific parameter array
  priorneutral = Prior([Uniform(priormu...),
     Uniform(priorcm...),
     Uniform(priorcellularity...)])

  priort = [2.0, log(Nmax)/log(2) + (MathConstants.eulergamma / log(2))]
  maxsel = selection(log(2), 0.99, priort[2], priort[2] - 1)
  priorsel = [0.0, maxsel]

priorselection = Prior([Uniform(priormu...),
   Uniform(priorcm...),
   Uniform(priorsel...),
   Uniform(priort...),
   Uniform(priorcellularity...)])

priorselection2 = Prior([Uniform(priormu...),
   Uniform(priorcm...),
   Uniform(priorsel...),
   Uniform(priort...),
   Uniform(priorsel...),
   Uniform(priort...),
   Uniform(priorcellularity...)])

  nparamatersneutral = 3
  nparametersselection = 5
  nparametersselection2 = 7

  #for the ABCSMC Model selection, each input is an array of model specific parameters,

  if maxclones == 1
    setup = ABCSMCModel(
    [tumourABCneutral, tumourABCselection], #simulation functions
    [nparamatersneutral, nparametersselection], # number of parameters for each model
    ϵT, #desired final tolerence
    [priorneutral, priorselection], #priors
    constants = [cst, cst], #constants
    nparticles = nparticles,
    maxiterations = maxiterations,
    convergence = convergence,
    modelkern = modelkern,
    ϵ1 = ϵ1,
    other = Nmaxinf);
  elseif maxclones == 2
    setup = ABCSMCModel(
    [tumourABCneutral, tumourABCselection, tumourABCselection2], #simulation functions
    [nparamatersneutral, nparametersselection, nparametersselection2], # number of parameters for each model
    ϵT, #desired final tolerence
    [priorneutral, priorselection, priorselection2], #priors
    constants = [cst, cst, cst], #constants
    nparticles = nparticles,
    maxiterations = maxiterations,
    convergence = convergence,
    modelkern = modelkern,
    ϵ1 = ϵ1,
    other = Nmaxinf);
  end

  return setup
end

"""
    fitABCmodels(data::Array{Float64, 1}, sname::String; <keyword arguments>)

Fit a stochastic model of cancer evolution to cancer sequencing data using Approximate Bayesian computation and infer the number of subclones (up to 2) and their relative fitness and time they emerge.
...
## Arguments
- `read_depth = 200.0`: Mean read depth of the target data set
- `minreads = 5`: Minimum number of reads needed to identify a mutation, this should be the same as was set in the sequencing data analysis
- `minvaf = 0.0`: Minimum VAF to identify a mutation, will override minreads if > 0.0
- `fmin = 0.01`: Minimum range of VAF to perform inference, if minvaf > 0.0 then fmin = minvaf
- `fmax = 0.75`: Maximum range of VAF to perform inference
- `maxiterations = 10^6`: Maximum number of iterations before ABC terminates
- `maxclones = 2`: Maximum number of clones, can be 1 or 2
- `nparticles = 500`: Number of particles (ie samples) in the ABC output
- `Nmax = 10^4`: Maximum population size used to fit data, increase if suspect that there is a late arising subclone that dominates the distribution
- `resultsdirectory = "output" `: Directory where posterior distributions will be saved to, this directory will be created if it does not already exist.
- `progress = true`: Show progress of ABC sampler with `ProgressMeter` package
- `verbose = true`: Print out summary at each ABC population step.
- `save = false `: Save output or not
- `ϵ1 = 10^6 `: Target ϵ for first ABC step
- `firstpass = false`: If set to true will run a limited first pass of the algorithm to determine a good starting ϵ1 if this unknown, otherwise ϵ1 is effectively infinity.
- `Nmaxinf = 10^10`: Scales selection coefficient value assuming the tumour size is Nmaxinf. Once value >10^9 has limited effect.
- `ρ = 0.0`: Overdispersion parameter for beta-binomial model of sequencing data. ρ = 0.0 means model is binomial sampling
- `adaptpriors = true`: If true priors on μ and clonalmutations are adapted based on the number of mutations in the data set which provide an upper an lower limit, this is an experimental feature that needs further validation, although initial tests suggest it performs well and does not skew inferences
- `timefunction = timefunc`: Function for KMC algorithm timestep. timefunc returns 1 meaning the timestep is the average of stochastic process. Alternatively timefuncrand can be specified which uses `-log(rand())` to increase the time step, so it is exponentially distributed rather than the mean of the exponential distribution. We use population doublings as our unit of time so this does not change the method and is an algorithmic choice.
- `ploidy = 2`: ploidy of the genome
- `d = 0.0`: Death rate of the the host population in the tumour, it is the compound parameter μ/β (where β = (b-d)/b) that can be inferred from the VAF distribution hence setting d = 0.0 means there is one fewer parameter in the ABC and hence increases efficiency.
- `b = log(2)`: Birth rate of the population. Default is set to `log(2)` so that tumour doubles with each unit increase in t in the absence of cell death.
- `mincellularity = 0.1`: If some prior knowledge on cellularity is known this can be modifed
- `maxcellularity = 1.1`: If some prior knowledge on cellularity is known this can be modifed. This is set to > 1.0 so that there is some flexibility in the inference.
- `convergence = 0.07`: Convergence for ABC. If new population ϵ is within convergence % of previous population then ABC stops.
- `savepopulations = false`: Save results from intermediate populations.
...
"""
function fitABCmodels(data::Array{Float64, 1}, sname::String;
  read_depth = 200.0, minreads = 5, minvaf = 0.0, fmin = 0.05,
  fmax = 0.75, maxiterations = 10^4, maxclones = 2,
  nparticles = 500, Nmax = 10^4, resultsdirectory::String = "output",
  progress = true, verbose = false, save = false,
  ϵ1 = 10^6, mincellularity = 0.1, maxcellularity = 1.1, firstpass = false,
  Nmaxinf = 10^10, ρ = 0.0,
  adaptpriors = true, timefunction = timefunc, ploidy = 2, d = 0.0, b = log(2),
  maxmu = 500, maxclonalmutations = 5000, convergence = 0.07,
  savepopulations = false)
  
  println("It is recommended that you use mobster for this type of analysis moving forward. mobster is an R package that provides similar functionality with orders of magnitude increases in speed and has many other features. SubClonalSelection.jl will remain here but I am unlikely to actively develop the package or be able to provide much ongoing support.")

  #make output directories
  if save
    makedirectories(joinpath(resultsdirectory, sname, "finalpopulation"))
  end
  #sequencing error is 1%, otherwise detection limit is controlled by minimum number of reads to call a variant
  detectionlimit = maximum([0.01, minreads/read_depth])

  #if a hard cutoff for minvaf is specified detectionlimit is changed also changes fmin in this case as its assumed all mutations < minvaf have been removed or shouldn't be used
  if minvaf > 0.0
    detectionlimit = minvaf
  end

  targetdata, VAF = gettargetDF(data, fmin = fmin, fmax = fmax)
  targetdataCDF = targetdata
  if save != false
    writedlm(joinpath(joinpath(resultsdirectory, sname, "finalpopulation"),
    "data", "$(sname).txt"), VAF)
  end

  dl = detectionlimit
  eps1 = ϵ1

  #adaptively set priors on μ and clonalmutations based on number of mutations
  # two extreme cases are mutations are all subclonal which defines a maximum mu,
  #or mutations are all clonal hich defines a maximum clonalmutational load
  if adaptpriors == true
    maxmu = 1.5 * (length(data) / ((1/detectionlimit) - 1))
    maxclonalmutations = length(data)
  end

  nparts = nparticles
  if (firstpass == true) && (length(VAF) < 3000) #I found this not to be useful when number of mutations is large
    println("################################################")
    println("Running first pass to get starting point for ABC")
    println("################################################")
    println("")
    ABCsetup = getsetup(1, detectionlimit = detectionlimit,
    read_depth = read_depth,
    nparticles = 100,
    maxiterations = 3 * nparts,
    modelkern = 0.5,
    Nmax = Nmax,
    mincellularity = mincellularity,
    maxcellularity = maxcellularity,
    maxclonalmutations = maxclonalmutations,
    maxmu = maxmu,
    timefunction = timefunction,
    ploidy = ploidy,
    d = d,
    b = b,
    fmin = fmin,
    fmax = fmax,
    Nmaxinf = Nmaxinf
    )
    abcres = runabcCancer(ABCsetup, targetdataCDF, verbose = verbose, progress = progress);
    eps1 = abcres.ϵ[maximum([1, length(abcres.ϵ) - 1])]
    println("################################################")
    println("Now running inference with ϵ1 = $(eps1)")
    println("################################################")
    println("")
  end

  ABCsetup = getsetup(maxclones, detectionlimit = dl,
  read_depth = read_depth,
  maxiterations = maxiterations,
  nparticles = nparticles,
  modelkern = 0.5,
  convergence = convergence,
  ϵ1 = eps1,
  Nmax = Nmax,
  mincellularity = mincellularity,
  ρ = ρ,
  maxclonalmutations = maxclonalmutations,
  maxmu = maxmu,
  timefunction = timefunction,
  ploidy = ploidy,
  d = d,
  b = b,
  fmin = fmin,
  fmax = fmax,
  Nmaxinf = Nmaxinf
  )

  abcres = runabcCancer(ABCsetup, targetdataCDF, verbose = verbose,
  progress = progress, Nmaxinf = Nmaxinf,
  savepopulations = savepopulations, sname = sname, resultsdirectory = resultsdirectory, VAF = VAF);

  posteriors, DFmp = getresults(abcres, joinpath(resultsdirectory), sname, VAF, save = save, Nmaxinf = Nmaxinf, savepopulations = false);

  finalresults = Results(ABCsetup, abcres, VAF, posteriors, DFmp, sname);

  return finalresults
end

"""
    fitABCmodels(data::String, sname::String; <keyword arguments>)

If data is a string will read in file. File should be a 1 column text file with VAF values in the rows.
"""
function fitABCmodels(data::String, sname::String;
  read_depth = 200.0, minreads = 5, fmin = 0.05, minvaf = 0.0,
  fmax = 0.75, maxiterations = 10^4, maxclones = 2,
  nparticles = 500, Nmax = 10^4, resultsdirectory::String = "output",
  progress = true, verbose = false, save = false,
  ϵ1 = 10^6, mincellularity = 0.1, maxcellularity = 1.1,
  firstpass = false, Nmaxinf = 10^10, ρ = 0.0,
  adaptpriors = true, timefunction = timefunc, ploidy = 2, d = 0.0, b = log(2),
  convergence = 0.07,
  savepopulations = false)

  VAF = readdlm(data)[:, 1]

  return fitABCmodels(VAF, sname;
  fmin = fmin, fmax = fmax, minreads = minreads, minvaf = minvaf,
  read_depth = read_depth, maxiterations = maxiterations,
  maxclones = maxclones, nparticles = nparticles,
  Nmax = Nmax, resultsdirectory = resultsdirectory,
  progress = progress, verbose = verbose, save = save,
  ϵ1 = ϵ1, mincellularity = mincellularity, maxcellularity = maxcellularity,
  firstpass = firstpass, Nmaxinf = Nmaxinf,
  ρ = ρ, adaptpriors = adaptpriors,
  timefunction = timefunction, ploidy = ploidy, d = d, b = b,
  convergence = convergence,
  savepopulations = savepopulations)
end
