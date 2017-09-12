function gettargetDF(VAF; fmin = 0.05, fmax = 0.75)

  targetdataDF = cumulativedistribution(VAF, fmin = fmin, fmax = fmax)

  return targetdataDF, VAF
end

function selection(λ, f, tend, t1)

    s = (λ .* t1 + log(f ./ (1 - f))) ./ (λ .* (tend - t1))
    return s
end


function getmodel(abcres, model)

    indeces = map(x -> x.model, abcres.particles) .== model
    abcresnew = deepcopy(abcres)
    abcresnew.particles = abcresnew.particles[indeces]
    abcresnew.parameters = abcresnew.parameters[model]

    return abcresnew
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
                Nmax = cst[6],
                timefunction = cst[7],
                detectionlimit = cst[8],
                s = Float64[],
                tevent = Float64[],
                cellularity = parameters[3])

    AD = CancerSeqSim.cumulativedist(simdata,
                    fmin = minimum(targetdata[:v]),
                    fmax = maximum(targetdata[:v]))

    return euclidean(AD.DF[:cumsum], targetdata[:cumsum]), [simdata.sampleddata.DF], true
end


function tumourABCselection(parameters, constants, targetdata)

    # we only consider

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

    AD = CancerSeqSim.cumulativedist(simdata,
                        fmin = minimum(targetdata[:v]),
                        fmax = maximum(targetdata[:v]))

    if length(simdata.output.subclonemutations) == 0
        out = [simdata.sampleddata.DF, 0.0, 0.0, 0.0]
    else
        out = [simdata.sampleddata.DF, simdata.output.subclonemutations[1], simdata.output.Ndivisions[1], simdata.output.clonefreq[1]]
    end

    c = false
    if ((sum(simdata.output.clonefreq.<0.95).==1) & (sum(simdata.output.clonefreq.>0.05).==1))[1] == true
        c = true
    end

    euclidean(AD.DF[:cumsum], targetdata[:cumsum]), out, c
end

function tumourABCselection2(parameters, constants, targetdata)

    # we only consider

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

    AD = CancerSeqSim.cumulativedist(simdata,
                        fmin = minimum(targetdata[:v]),
                        fmax = maximum(targetdata[:v]))

    if length(simdata.output.subclonemutations) < 2
        out = [simdata.sampleddata.DF, 0.0, 0.0, 0.0]
    else
        out = [simdata.sampleddata.DF, simdata.output.subclonemutations[1],
        simdata.output.subclonemutations[2],
        simdata.output.Ndivisions[1],
        simdata.output.Ndivisions[2],
        simdata.output.clonefreq[1],
        simdata.output.clonefreq[2]]
    end

    c = false
    if ((sum(simdata.output.clonefreq.<0.95).==2) & (sum(simdata.output.clonefreq.>0.05).==2))[1] == true
        c = true
    end

    euclidean(AD.DF[:cumsum], targetdata[:cumsum]), out, c
end

function getsetup(maxclones; nparticles = 100, maxiterations = 10^4, convergence = 0.1, ϵT = 1.0, read_depth = 100.0, Nmax = 10^3, detectionlimit = 0.05, modelkern = 0.4, scalefactor = 6, ϵ1 = 10^6, mincellularity = 0.1, ρ = 0.0)

  #function to define priors, constants and create ABCSMC model type

  timefunc() = 1
  ploidy = 2
  read_depth = read_depth
  d = 0.0
  b = log(2)
  Nmax = Nmax
  ρ = ρ

  cst = [ploidy, read_depth, d, b, ρ, Nmax, timefunc, detectionlimit];

  #priors
  #max mu is 1e-7 per bp per division, max cm is
  priormu = [0.01, 620.0]
  priorcm = [0.0, 31000.0]
  priorcellularity = [mincellularity, 1.0]

  #need to create Prior type which has a distribution type array with a corresponding distribution specific parameter array
  priorneutral = Prior([Uniform(priormu...),
     Uniform(priorcm...),
     Uniform(priorcellularity...)])

  priorsel = [3.0, 25.0]
  priort = [3.0, log(Nmax)/log(2) + (eulergamma / log(2))]

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
    scalefactor = scalefactor,
    ϵ1 = ϵ1);
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
    scalefactor = scalefactor,
    ϵ1 = ϵ1);
  end

  return setup
end

"""
    fitABCmodels(data::Array{Float64, 1}, sname::String; <keyword arguments>)

Fit a stochastic model of cancer evolution to cancer sequencing data using Approximate Bayesian computation and infer the number of subclones (up to 2) and their relative fitness and time they emerge.
...
## Arguments
- `read_depth = 200.0`: Mean read depth of the target data set
- `minreads = 5`: Minimum number of reads to identify a mutation
- `fmin = 0.01`: Minimum range of VAF to perform inference
- `fmax = 0.75`: Maximum range of VAF to perform inference
- `maxiterations = 10^4 `: Maximum number of iterations before ABC terminates
- `maxclones = 2`: Maximum number of clones, can be 1 or 2
- `nparticles = 500`: Number of particles (ie samples) in the ABC output
- `Nmax = 10^4`: Maximum population size used to fit data, increase if suspect that there is a late arising clone
- `resultsdirectory = "output" `: Directory where posterior distributions will be saved to
- `progress = true`: Show progress of ABC sampler with `ProgressMeter` package
- `verbose = true`: Print out summary at each ABC population step.
- `save = false `: Save output or not
- `inferdetection = false `: Do a first past to infer cellularity to modify sequencing depth and detection limit
- `ϵ1 = 10^6 `: Target ϵ for first ABC step, if you find the model with 2 clones often dies out, decrease this value
- `firstpass = false`: If set to true will run a limited first pass of the algorithm to determine a good starting ϵ1 if this unkown.
- `Nmaxinf = 10^10`: Scales selection coefficient value assuming the tumour size is Nmaxinf. Once value >10^9 has limited effect.
- `scalefactor = 6`: Parameter for perturbation kernel for parameter values. Larger values means space will be explored more slowly but fewer particles will be perturbed outside prior range.
- `ρ = 0.0`: Overdispersion parameter for beta-binomial model of sequencing data. ρ = 0.0 means model is binomial sampling
...
"""
function fitABCmodels(data::Array{Float64, 1}, sname::String; read_depth = 200.0, minreads = 5, fmin = 0.01, fmax = 0.75, maxiterations = 10^4, maxclones = 2, nparticles = 500, Nmax = 10^4, resultsdirectory::String = "output", progress = true, verbose = true, save = false, inferdetection = false, ϵ1 = 10^6, mincellularity = 0.1, firstpass = false, Nmaxinf = 10^10, scalefactor = 6, ρ = 0.0)

  #make output directories
  if save != false
    makedirectory(resultsdirectory)
    makedirectories(joinpath(resultsdirectory, sname))
  end
  #sequencing error is 1%, otherwise detection limit is controlled by minimum number of reads to call a variant
  detectionlimit = maximum([0.01, minreads/read_depth])

  targetdata, VAF = gettargetDF(data, fmin = fmin, fmax = fmax)
  targetdataDF = targetdata.DF
  if save != false
    writedlm(joinpath(joinpath(resultsdirectory, sname), "data", "$(sname).txt"), VAF)
  end

  dl = detectionlimit
  eps1 = ϵ1

  if inferdetection == true
    abcsetup = getsetup(1, detectionlimit = detectionlimit,
    read_depth = read_depth,
    maxiterations = 10^3,
    nparticles = 100,
    modelkern = 0.5,
    scalefactor = 6,
    Nmax = Nmax,
    mincellularity = mincellularity
    )
    abcres = ApproxBayes.runabcCancer(abcsetup, targetdataDF, verbose = verbose, progress = progress);
    DFpost0 = collectoutput0clone(getmodel(abcres, 1))
    c = mean(DFpost0[:cellularity])

    println()
    println("####################################################")
    println("Now running real ABC with detection limit $(5./(c*read_depth))")
    println("####################################################")
    println()
    dl = 5./(c*read_depth)
    eps1 = abcres.ϵ[end]
  end

  nparts = nparticles
  if (firstpass == true) && (length(VAF) < 3000)
    println("################################################")
    println("Running first pass to get starting point for ABC")
    println("################################################")
    println("")
    abcsetup = getsetup(1, detectionlimit = detectionlimit,
    read_depth = read_depth,
    nparticles = 100,
    maxiterations = 3 * nparts,
    modelkern = 0.5,
    scalefactor = scalefactor,
    Nmax = Nmax,
    mincellularity = mincellularity
    )
    abcres = ApproxBayes.runabcCancer(abcsetup, targetdataDF, verbose = verbose, progress = progress);
    eps1 = abcres.ϵ[maximum([1, length(abcres.ϵ) - 1])]
    println("################################################")
    println("Now running inference with ϵ1 = $(eps1)")
    println("################################################")
    println("")
  end

  abcsetup = getsetup(maxclones, detectionlimit = dl,
  read_depth = read_depth,
  maxiterations = maxiterations,
  nparticles = nparticles,
  modelkern = 0.5,
  scalefactor = scalefactor,
  convergence = 0.08,
  ϵ1 = eps1,
  Nmax = Nmax,
  mincellularity = mincellularity,
  ρ = ρ
  )
  abcres = ApproxBayes.runabcCancer(abcsetup, targetdataDF, verbose = verbose, progress = progress);

  show(abcres)

  posteriors, DFmp = getresults(abcres, joinpath(resultsdirectory, sname), sname, VAF, save = save, Nmaxinf = Nmaxinf)

  return Results(abcsetup, abcres, VAF, posteriors, DFmp, sname)
end

"""
    fitABCmodels(data::String, sname::String; <keyword arguments>)

If data is a string will read in file. File should be a 1 column text file with VAF values in the rows.
"""
function fitABCmodels(data::String, sname::String; read_depth = 200.0, minreads = 5, fmin = 0.01, fmax = 0.75, maxiterations = 10^4, maxclones = 2, nparticles = 500, Nmax = 10^3, resultsdirectory::String = "output", progress = true, verbose = true, save = false, inferdetection = false, ϵ1 = 10^6, mincellularity = 0.1, firstpass = true, Nmaxinf = 10^10, scalefactor = 6, ρ = 0.0)

  VAF = readdlm(data)[:, 1]

  return fitABCmodels(VAF, sname; fmin = fmin, fmax = fmax, minreads = minreads, read_depth = read_depth, maxiterations = maxiterations, maxclones = maxclones, nparticles = nparticles, Nmax = Nmax, resultsdirectory = resultsdirectory, progress = progress, verbose = verbose, save = save, inferdetection = inferdetection, ϵ1 = ϵ1, mincellularity = mincellularity, firstpass = firstpass, Nmaxinf = Nmaxinf, scalefactor = scalefactor, ρ = ρ)
end
