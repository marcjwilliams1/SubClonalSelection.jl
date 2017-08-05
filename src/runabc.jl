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

function getsetup(maxclones; nparticles = 100, maxiterations = 10^4, convergence = 0.1, ϵT = 1.0, read_depth = 100.0, Nmax = 10^3, detectionlimit = 0.05, modelkern = 0.5, scalefactor = 4, ϵ1 = 10^6)

  #function to define priors, constants and create ABCSMC model type

  timefunc() = 1
  ploidy = 2
  read_depth = read_depth
  d = 0.0
  b = log(2)
  Nmax = Nmax
  ρ = 0.0

  cst = [ploidy, read_depth, d, b, ρ, Nmax, timefunc, detectionlimit];

  #priors
  priormu = [0.01, 200.0]
  priorcm = [0.0, 5000.0]
  priorcellularity = [0.1, 1.0]

  #need to create Prior type which has a distribution type array with a corresponding distribution specific parameter array
  priorneutral = Prior([Uniform(priormu...),
     Uniform(priorcm...),
     Uniform(priorcellularity...)])

  priorsel = [0.0, 20.0]
  priort = [3.0, 12.0]

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

function fitABCmodels(data::Array{Float64, 1}, sname; fmin = 0.05, fmax = 0.75, detectionlimit = fmin, read_depth = 200.0, maxiterations = 10^4, maxclones = 2, nparticles = 500, Nmax = 10^3, resultsdirectory::String = "output", progress = true, verbose = true, save = false, inferdetection = false, ϵ1 = 10^6)

  #make output directories
  if save != false
    makedirectories(joinpath(resultsdirectory, sname))
  end

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
    scalefactor = 4,
    Nmax = Nmax
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

  abcsetup = getsetup(maxclones, detectionlimit = dl,
  read_depth = read_depth,
  maxiterations = maxiterations,
  nparticles = nparticles,
  modelkern = 0.5,
  scalefactor = 4,
  convergence = 0.05,
  ϵ1 = eps1,
  Nmax = Nmax
  )
  abcres = ApproxBayes.runabcCancer(abcsetup, targetdataDF, verbose = verbose, progress = progress);

  show(abcres)

  posteriors, DFmp = getresults(abcres, joinpath(resultsdirectory, sname), sname, save = save)

  return Results(abcsetup, abcres, VAF, posteriors, DFmp)
end

function fitABCmodels(data::String, sname; fmin = 0.05, fmax = 0.75, detectionlimit = fmin, read_depth = 200.0, maxiterations = 10^4, maxclones = 2, nparticles = 500, Nmax = 10^3, resultsdirectory::String = "output", progress = true, verbose = true, save = false, inferdetection = false, ϵ1 = 10^6)

  VAF = readdlm(data)[:, 1]

  return fitABCmodels(VAF, sname; fmin = fmin, fmax = fmax, detectionlimit = detectionlimit, read_depth = read_depth, maxiterations = maxiterations, maxclones = maxclones, nparticles = nparticles, Nmax = Nmax, resultsdirectory = resultsdirectory, progress = progress, verbose = verbose, save = save, inferdetection = inferdetection, ϵ1 = ϵ1)
end
