Mcdf(f,fmin,fmax) = (1.0./f - 1.0/fmax) ./ (1.0/fmin - 1.0/fmax)

function selection(λ, f, tend, t1)
    #define the equation for selection as above
    s = (λ .* t1 + log.(f ./ (1 - f))) ./ (λ .* (tend - t1))
    return s
end

function selection2clone(λ, f1, f2, tend, t1, t2)
    #define the equation for selection as above

    s1 = zeros(Float64, length(f1))
    s2 = zeros(Float64, length(f1))

    for i in 1:length(f1)
      if (f2[i] + f1[i]) < 1.0
        s1[i] = (λ .* t1[i] + log.(f1[i] ./ (1 - f1[i] - f2[i]))) ./ (λ .* (tend[i] - t1[i]))
        s2[i] = (λ .* t2[i] + log.(f2[i] ./ (1 - f1[i] - f2[i]))) ./ (λ .* (tend[i] - t2[i]))
      else
        s1[i] = (λ .* t1[i] + log.((f1[i] - f2[i]) ./ (1 - f1[i]))) ./ (λ .* (tend[i] - t1[i]))
        s2[i] = (λ .* t2[i] + log.(f2[i] ./ (1 - f1[i]))) ./ (λ .* (tend[i] - t2[i]))
      end
    end

    return s1, s2
end

function collectoutput1clone(abcres; Nmax = 10^10)

    scmuts = map(x -> x.other[2], abcres.particles)
    scdivs = map(x -> x.other[3], abcres.particles)
    scfreq = map(x -> x.other[4], abcres.particles)
    mu = abcres.parameters[:, 1]
    weights = abcres.weights

    # eulergamma/log(2) is stochastic correction see Durrett Branching Process Models of Cancer, needed for selection calculation
    t1 = ((shuffle(scmuts) ./ mu) / (2 * log.(2))) - eulergamma/log.(2) #shuffle so that t1 posterior is not correlated
    tend = (log.(Nmax .* (1 - scfreq)) / log.(2))
    s = selection(log.(2), scfreq, tend, t1)

    DF = DataFrame(mu = mu,
    clonalmutations = abcres.parameters[:, 2],
    s = s,
    t = t1,
    cellularity = abcres.parameters[:, 5],
    frequency = scfreq,
    scmuts = map(x -> Float64(x), scmuts),
    weight = weights)

    return DF
end

function swapvalues(x1, x2, indeces)

  swap1 = x1[indeces]
  swap2 = x2[indeces]
  x1[indeces] = swap2
  x2[indeces] = swap1

  return x1, x2
end

function clonesize(scfreq1, scfreq2)
  #assume that if scfreq1 + scfreq2 > 1 then clones are nested
  freqfactor = zeros(Float64, length(scfreq1))
  for i in 1:length(scfreq1)
    if (scfreq1[i] + scfreq2[i]) .> 1.0
      freqfactor[i] = 1 - scfreq1[i]
      freqfactor[i] = maximum([0.001, freqfactor[i]])
    else
      freqfactor[i] = 1 - scfreq1[i] - scfreq2[i]
      freqfactor[i] = maximum([0.001, freqfactor[i]])
    end
  end

  return freqfactor
end

function collectoutput2clone(abcres; Nmax = 10^10)

    scmuts1 = map(x -> x.other[2], abcres.particles)
    scmuts2 = map(x -> x.other[3], abcres.particles)
    scdivs1 = map(x -> x.other[4], abcres.particles)
    scdivs2 = map(x -> x.other[5], abcres.particles)
    scfreq1 = map(x -> x.other[6], abcres.particles)
    scfreq2 = map(x -> x.other[7], abcres.particles)
    mu = abcres.parameters[:, 1]

    s1 = abcres.parameters[:, 3]
    t1 = abcres.parameters[:, 4]
    s2 = abcres.parameters[:, 5]
    t2 = abcres.parameters[:, 6]
    weights = abcres.weights

    scmuts1 = map(x -> Float64(x), scmuts1)
    scmuts2 = map(x -> Float64(x), scmuts2)

    indeces = .!(scfreq1 .> scfreq2)

    #identify subclone 1 as the highest frequency subclone
    scfreq1, scfreq2 = swapvalues(scfreq1, scfreq2, indeces)
    scmuts1, scmuts2 = swapvalues(scmuts1, scmuts2, indeces)
    scdivs1, scdivs2 = swapvalues(scdivs1, scdivs2, indeces)
    s1, s2 = swapvalues(s1, s2, indeces)
    t1, t2 = swapvalues(t1, t2, indeces)

    freqfactor = clonesize(scfreq1, scfreq2)

    t1a = ((shuffle(scmuts1) ./ mu) / (2 * log.(2))) - eulergamma/log.(2)
    t1b = ((shuffle(scmuts2) ./ mu) / (2 * log.(2))) - eulergamma/log.(2)
    tend = (log.(Nmax .* (freqfactor)) / log.(2))
    s1, s2 = selection2clone(log.(2), scfreq1, scfreq2, tend, t1a, t1b)

    DF = DataFrame(mu = abcres.parameters[:, 1],
    clonalmutations = abcres.parameters[:, 2],
    s1 = s1,
    t1 = t1a,
    s2 = s2,
    t2 = t1b,
    cellularity = abcres.parameters[:, 7],
    frequency1 = scfreq1,
    frequency2 = scfreq2,
    scmuts1 = scmuts1,
    scmuts2 = scmuts2,
    weight = weights)

    return DF
end

function collectoutput0clone(abcres)

  mupost = abcres.parameters[:, 1]
  cmpost = abcres.parameters[:, 2]
  cellularity = abcres.parameters[:, 3]
  weights = abcres.weights

  DFpost = DataFrame(mu = mupost,
  clonalmutations = cmpost,
  cellularity = cellularity,
  weight = weights)
end

function cumulativedistribution(VAF; fmin = 0.1, fmax = 0.3)

    #calculate cumulative sum
    steps = fmax:-0.001:fmin
    cumsum = Array{Int64}(0)
    v = Array{Float64}(0)

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

function averagehistogram(particles, model, VAF)

    x = 0.005:0.01:1.005
    y = fit(Histogram, VAF, x, closed=:right)
    DFhist = DataFrame(VAF = x[1:end-1], freq = y.weights)

    particles = particles[map(x -> x.model, particles).==model]
    wts = map(x -> x.weight, particles)
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

    for i in 1:size(M, 1)
      push!(mvalues, median(vec(collect(M[i, :]')), weights(wts)))
      push!(meanvalues, mean(vec(collect(M[i, :]')), weights(wts)))
      push!(lquant, quantile(vec(collect(M[i, :]')), weights(wts), 0.025))
      push!(uquant, quantile(vec(collect(M[i, :]')), weights(wts), 0.975))
      push!(lquart, quantile(vec(collect(M[i, :]')), weights(wts), 0.25))
      push!(uquart, quantile(vec(collect(M[i, :]')), weights(wts), 0.75))
    end

    DFr = DataFrame(median = mvalues,
                    mean = meanvalues,
                    lowerq95 = lquant,
                    upperq95 = uquant,
                    lowerquartile = lquart,
                    upperquartile = uquart,
                    VAF = DFhist[:VAF],
                    truecounts = DFhist[:freq])

    return DFr
end

function saveresults(res::Results; resultsdirectory = "output")
  makedirectories(joinpath(resultsdirectory, res.SampleName))
  getresults(res.ABCresults, resultsdirectory, res.SampleName, res.VAF, save = true)
  return
end

function getresults(abcres, resultsdirectory, sname, VAF; save = false, Nmaxinf = 10^10)

  resultsdirectory = joinpath(resultsdirectory, sname)
  posteriors = Posterior[]
  #save model posterior
  DFmp = DataFrame(Model = map(x -> string(x),0:length(abcres.modelprob) - 1), Probability = abcres.modelprob)

  if abcres.modelprob[1] > 0.0
    DFpost0 = collectoutput0clone(getmodel(abcres, 1))
    DFr = averagehistogram(abcres.particles, 1, VAF)
    if save == true
      CSV.write(joinpath(resultsdirectory, "posterior", "$(sname)-parameters-clone0.csv"), DFpost0)
      CSV.write(joinpath(resultsdirectory, "posterior", "$(sname)-histogram-clone0.csv"), DFr)
    end
    push!(posteriors, Posterior(DFr, DFpost0, abcres.modelprob[1]))
  else
    push!(posteriors, Posterior(DataFrame(), DataFrame(), abcres.modelprob[1]))
  end

  if abcres.modelprob[2] > 0.0
    DFpost1 = collectoutput1clone(getmodel(abcres, 2), Nmax = Nmaxinf)
    DFr = averagehistogram(abcres.particles, 2, VAF)
    if save == true
      CSV.write(joinpath(resultsdirectory, "posterior", "$(sname)-histogram-clone1.csv"), DFr)
      CSV.write(joinpath(resultsdirectory, "posterior", "$(sname)-parameters-clone1.csv"), DFpost1)
    end
    push!(posteriors, Posterior(DFr, DFpost1, abcres.modelprob[2]))
  else
    push!(posteriors, Posterior(DataFrame(), DataFrame(), abcres.modelprob[2]))
  end

  if (length(abcres.modelprob) > 2) && (abcres.modelprob[3] > 0.0)
    DFpost2 = collectoutput2clone(getmodel(abcres, 3), Nmax = Nmaxinf)
    DFr = averagehistogram(abcres.particles, 3, VAF)
    if save == true
      CSV.write(joinpath(resultsdirectory, "posterior", "$(sname)-parameters-clone2.csv"), DFpost2)
      CSV.write(joinpath(resultsdirectory, "posterior", "$(sname)-histogram-clone2.csv"), DFr)
    end
    push!(posteriors, Posterior(DFr, DFpost2, abcres.modelprob[3]))
  else
    if (length(abcres.modelprob) > 2)
      push!(posteriors, Posterior(DataFrame(), DataFrame(), abcres.modelprob[3]))
    end
  end

  DF = DataFrame(Model = map(x -> string(x),0:length(abcres.modelprob) - 1), Probability = abcres.modelprob)
  if save == true
    CSV.write(joinpath(resultsdirectory, "posterior", "$(sname)-modelprobabilities.csv"), DF)
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

function makedirectory(resultsdirectory)
  if isdir(joinpath(resultsdirectory)) == false
    mkdir(joinpath(resultsdirectory))
  end
end

function show(io::IO, res::Results)
  show(res.ABCresults)
end

function show(io::IO, ABCresults::ApproxBayes.ABCSMCmodelresults)

  posteriors, DFmp = getresults(ABCresults, "nothing", "nothing", [1.0, 2.0]; save = false, Nmaxinf = ABCresults.setup.other)

  @printf("Total number of simulations: %.2e\n", sum(ABCresults.numsims))
  println("Cumulative number of simulations = $(cumsum(ABCresults.numsims))")
  @printf("Acceptance ratio: %.2e\n\n", ABCresults.accratio)
  println("Tolerance schedule = $(round.(ABCresults.ϵ, 2))\n")

  print("Model probabilities:\n")
  for j in 1:length(ABCresults.modelprob)
    @printf("\tModel %d (%d subclones): %.3f\n", j, j-1, ABCresults.modelprob[j])
  end

  print("\nParameters:\n\n")

  for j in 1:length(ABCresults.parameters)
    if ABCresults.modelprob[j] > 0.0
      print("Model $j ($(j-1) subclones)\n")

      upperci = zeros(Float64, size(posteriors[j].Parameters, 2) - 1)
      lowerci = zeros(Float64, size(posteriors[j].Parameters, 2) - 1)
      parametermeans = zeros(Float64, size(posteriors[j].Parameters, 2) - 1)
      parametermedians = zeros(Float64, size(posteriors[j].Parameters, 2) - 1)

      for i in 1:(size(posteriors[j].Parameters, 2) - 1)
        parametermeans[i] = mean(posteriors[j].Parameters[:, i],
        weights(posteriors[j].Parameters[:weight]))
        parametermedians[i] = median(posteriors[j].Parameters[:, i],
        weights(posteriors[j].Parameters[:weight]))
        (lowerci[i], upperci[i]) = quantile(posteriors[j].Parameters[:, i],
        weights(posteriors[j].Parameters[:weight]), [0.025,0.975])
      end

      print("\tMedian (95% intervals):\n")
      for i in 1:length(parametermeans)
          if j == 1
            parameternames = ["μ/β", "Clonal Mutations", "Cellularity"]
          elseif j == 2
            parameternames = ["μ/β", "Clonal Mutations", "Fitness", "Time (tumour doublings)", "Cellularity", "Subclone Frequency", "Subclone Mutations"]
          elseif j == 3
            parameternames = ["μ/β", "Clonal Mutations", "Fitness - Subclone 1", "Time (tumour doublings) - Subclone 1", "Fitness - Subclone 2", "Time (tumour doublings) - Subclone 2", "Cellularity", "Subclone 1 Frequency", "Subclone 2 Frequency", "Subclone 1 Mutations", "Subclone 2 Mutations"]
          end
          @printf("\tParameter %d - %s: %.2f (%.2f,%.2f)\n", i, parameternames[i], parametermedians[i], lowerci[i], upperci[i])
      end
    end
  end
end
