Mcdf(f,fmin,fmax) = (1.0./f - 1.0/fmax) ./ (1.0/fmin - 1.0/fmax)


function collectoutput1clone(abcres)

    scmuts = map(x -> x.other[2], abcres.particles)
    scdivs = map(x -> x.other[3], abcres.particles)
    scfreq = map(x -> x.other[4], abcres.particles)

    DF = DataFrame(mu = abcres.parameters[:, 1],
    clonalmutations = abcres.parameters[:, 2],
    s = abcres.parameters[:, 3],
    t = abcres.parameters[:, 4],
    cellularity = abcres.parameters[:, 5],
    freq = scfreq,
    scmuts = map(x -> Float64(x), scmuts))

    return DF
end

function swapvalues(x1, x2, indeces)

  swap1 = x1[indeces]
  swap2 = x2[indeces]
  x1[indeces] = swap2
  x2[indeces] = swap1

  return x1, x2
end

function collectoutput2clone(abcres)

    scmuts1 = map(x -> x.other[2], abcres.particles)
    scmuts2 = map(x -> x.other[3], abcres.particles)
    scdivs1 = map(x -> x.other[4], abcres.particles)
    scdivs2 = map(x -> x.other[5], abcres.particles)
    scfreq1 = map(x -> x.other[6], abcres.particles)
    scfreq2 = map(x -> x.other[7], abcres.particles)

    s1 = abcres.parameters[:, 3]
    t1 = abcres.parameters[:, 4]
    s2 = abcres.parameters[:, 5]
    t2 = abcres.parameters[:, 6]

    scmuts1 = map(x -> Float64(x), scmuts1)
    scmuts2 = map(x -> Float64(x), scmuts2)

    indeces = !(scfreq1 .> scfreq2)

    scfreq1, scfreq2 = swapvalues(scfreq1, scfreq2, indeces)
    scmuts1, scmuts2 = swapvalues(scmuts1, scmuts2, indeces)
    scdivs1, scdivs2 = swapvalues(scdivs1, scdivs2, indeces)
    s1, s2 = swapvalues(s1, s2, indeces)
    t1, t2 = swapvalues(t1, t2, indeces)

    DF = DataFrame(mu = abcres.parameters[:, 1],
    clonalmutations = abcres.parameters[:, 2],
    s1 = s1,
    t1 = t1,
    s2 = s2,
    t2 = t2,
    cellularity = abcres.parameters[:, 7],
    freq1 = scfreq1,
    freq2 = scfreq2,
    scmuts1 = scmuts1,
    scmuts2 = scmuts2)

    return DF
end

function collectoutput0clone(abcres)

  mupost = abcres.parameters[:, 1]
  cmpost = abcres.parameters[:, 2]
  cellularity = abcres.parameters[:, 3]

  DFpost = DataFrame(mu = mupost, clonalmutations = cmpost, cellularity = cellularity)
end

function collectoutput(abcres, tend)

    Nend = exp(tend*log(2))

    scmuts = map(x -> x.other[2], abcres.particles)
    scdivs = map(x -> x.other[3], abcres.particles)
    scfreq = map(x -> x.other[4], abcres.particles)

    DF = DataFrame(mu = abcres.parameters[:, 1],
    clonalmutations = abcres.parameters[:, 2],
    freq = scfreq,
    time = (1/(2*log(2))) .* (scmuts./abcres.parameters[:, 1] ),
    ndivs = scmuts./abcres.parameters[:, 1],
    scmuts = map(x -> Float64(x), scmuts))

    #euler gamma is a stochastic correction
    DF[:s] = selection(log(2), scfreq, (log((1 - scfreq) * Nend)/log(2) ) + + (eulergamma/log(2)), Array(DF[:time]))
    return DF
end

function cumulativedistribution(VAF; fmin = 0.1, fmax = 0.3)

    #calculate cumulative sum
    steps = fmax:-0.001:fmin
    cumsum = Array(Int64, 0)
    v = Array(Float64, 0)

    for i in steps
        push!(cumsum, sum(VAF .>= i))
        push!(v, i)
    end
    cumsum = cumsum - cumsum[1]

    DF = DataFrame(cumsum = map(Float64, cumsum), v = v)
    DF[:invf] = 1 ./ DF[:v] - 1 ./ fmax

    DF[:theory] = Mcdf(DF[:v], fmin, fmax)
    DF[:normalized] = DF[:cumsum] ./ maximum(DF[:cumsum])

    lmfit = fit(LinearModel, @formula(cumsum ~ invf + 0), DF)
    DF[:prediction] = predict(lmfit)

    return CancerSeqSim.AnalysedData(DF, VAF)
end


function averagehistogram(particles)

    N = length(particles)

    M = zeros(Int64, 100, N)
    i = 1

    for j in 1:N

        M[:, i] = freqcounts = convert(Array, particles[i].other[1][:freq])
        i = i + 1

    end

    mvalues = Float64[]
    meanvalues = Float64[]
    lquant = Float64[]
    lquart = Float64[]
    uquant = Float64[]
    uquart = Float64[]
    sd = Float64[]

    for i in 1:size(M, 1)
      push!(mvalues, median(vec(collect(M[i, :]'))))
      push!(meanvalues, mean(vec(collect(M[i, :]'))))
      push!(lquant, quantile(vec(collect(M[i, :]')), 0.025))
      push!(uquant, quantile(vec(collect(M[i, :]')), 0.975))
      push!(lquart, quantile(vec(collect(M[i, :]')), 0.25))
      push!(uquart, quantile(vec(collect(M[i, :]')), 0.75))
      push!(sd, std(vec(collect(M[i, :]'))))
    end

    DFr = DataFrame(median = mvalues,
                    mean = meanvalues,
                    lowerq = lquant,
                    upperq = uquant,
                    lowerquartile = lquart,
                    upperquartile = uquart,
                    sd = sd)

    return DFr
end

function averagehistogram(particles, model, VAF)

    x = 0.005:0.01:1.005
    y = fit(Histogram, VAF, x, closed=:right)
    DFhist = DataFrame(VAF = x[1:end-1], freq = y.weights)

    particles = particles[map(x -> x.model, particles).==model]
    N = length(particles)

    M = zeros(Int64, 100, N)
    i = 1

    for j in 1:N

        M[:, i] = convert(Array, particles[i].other[1][:freq])
        i = i + 1

    end

    mvalues = Float64[]
    meanvalues = Float64[]
    lquant = Float64[]
    lquart = Float64[]
    uquant = Float64[]
    uquart = Float64[]
    sd = Float64[]

    for i in 1:size(M, 1)
      push!(mvalues, median(vec(collect(M[i, :]'))))
      push!(meanvalues, mean(vec(collect(M[i, :]'))))
      push!(lquant, quantile(vec(collect(M[i, :]')), 0.025))
      push!(uquant, quantile(vec(collect(M[i, :]')), 0.975))
      push!(lquart, quantile(vec(collect(M[i, :]')), 0.25))
      push!(uquart, quantile(vec(collect(M[i, :]')), 0.75))
      push!(sd, std(vec(collect(M[i, :]'))))
    end

    DFr = DataFrame(median = mvalues,
                    mean = meanvalues,
                    lowerq = lquant,
                    upperq = uquant,
                    lowerquartile = lquart,
                    upperquartile = uquart,
                    sd = sd,
                    VAF = DFhist[:VAF],
                    truecounts = DFhist[:freq])

    return DFr
end

function saveresults(res::Results; resultsdirectory = "output")

  makedirectories(joinpath(resultsdirectory, sname))
  getresults(res.ABCresults, resultsdirectory, sname, save = true)
  return
end

function getresults(abcres, resultsdirectory, sname, VAF; save = false)

  posteriors = Posterior[]

  #save model posterior
  DFmp = DataFrame(Model = map(x -> string(x),0:length(abcres.modelprob) - 1), Probability = abcres.modelprob)

  if abcres.modelprob[1] > 0.0
    DFpost0 = collectoutput0clone(getmodel(abcres, 1))
    DFr = averagehistogram(abcres.particles, 1, VAF)
    if save == true
      writetable(joinpath(resultsdirectory, "posterior", "$(sname)-parameters-clone0.csv"), DFpost0)
      writetable(joinpath(resultsdirectory, "posterior", "$(sname)-histogram-clone0.csv"), DFr)
    end
    push!(posteriors, Posterior(DFr, DFpost0, abcres.modelprob[1]))
  else
    push!(posteriors, Posterior(DataFrame(), DataFrame(), abcres.modelprob[1]))
  end

  if abcres.modelprob[2] > 0.0
    DFpost1 = collectoutput1clone(getmodel(abcres, 2))
    DFr = averagehistogram(abcres.particles, 2, VAF)
    if save == true
      writetable(joinpath(resultsdirectory, "posterior", "$(sname)-histogram-clone1.csv"), DFr)
      writetable(joinpath(resultsdirectory, "posterior", "$(sname)-parameters-clone1.csv"), DFpost1)
    end
    push!(posteriors, Posterior(DFr, DFpost1, abcres.modelprob[2]))
  else
    push!(posteriors, Posterior(DataFrame(), DataFrame(), abcres.modelprob[2]))
  end

  if (length(abcres.modelprob) > 2) & (abcres.modelprob[3] > 0.0)
    DFpost2 = collectoutput2clone(getmodel(abcres, 3))
    DFr = averagehistogram(abcres.particles, 3, VAF)
    if save == true
      writetable(joinpath(resultsdirectory, "posterior", "$(sname)-parameters-clone2.csv"), DFpost2)
      writetable(joinpath(resultsdirectory, "posterior", "$(sname)-histogram-clone2.csv"), DFr)
    end
    push!(posteriors, Posterior(DFr, DFpost2, abcres.modelprob[3]))
  else
    push!(posteriors, Posterior(DataFrame(), DataFrame(), abcres.modelprob[3]))
  end

  DF = DataFrame(Model = map(x -> string(x),0:length(abcres.modelprob) - 1), Probability = abcres.modelprob)
  if save == true
    writetable(joinpath(resultsdirectory, "posterior", "$(sname)-modelprobabilities.csv"), DF)
  end

  return posteriors, DFmp
end

function makedirectories(resultsdirectory)

  if isdir(resultsdirectory) == false
    mkdir(resultsdirectory)
  end

  if isdir(joinpath(resultsdirectory, "plots")) == false
    mkdir(joinpath(resultsdirectory, "plots"))
  end
  if isdir(joinpath(resultsdirectory, "processed")) == false
    mkdir(joinpath(resultsdirectory, "processed"))
  end

  if isdir(joinpath(resultsdirectory, "posterior")) == false
    mkdir(joinpath(resultsdirectory, "posterior"))
  end


    if isdir(joinpath(resultsdirectory, "data")) == false
      mkdir(joinpath(resultsdirectory, "data"))
    end

end

function makeplotsdirectories(resultsdirectory)
  if isdir(joinpath(resultsdirectory)) == false
    mkdir(joinpath(resultsdirectory))
  end
  if isdir(joinpath(resultsdirectory, "plots")) == false
    mkdir(joinpath(resultsdirectory, "plots"))
  end
end

function show(res::Results)

  show(res.ABCresults)
end
