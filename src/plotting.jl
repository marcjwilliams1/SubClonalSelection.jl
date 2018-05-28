function plothistogram(res, model = 0; removelowfrequencies = true, plotsbackend = nothing)
  if plotsbackend != nothing
      plotsbackend()
  end
  if removelowfrequencies == true
      fmin = res.ABCsetup.Models[1].constants[9]
  else
      fmin = 0.0
  end

  model = model + 1

  DFres = res.Posterior[model].MeanHistogram
  DFres[:VAF] = collect(0.01:0.01:1.0)
  DF = DataFrame(VAF = res.VAF)
  DFres = DFres[1:75, :]
  DFreshist = deepcopy(DFres)
  DFres = DFres[DFres[:VAF].>=fmin, :]

  if model - 1 == 0
      postcolor = RGBA(0.0, 76/255, 153/255) # plot neutral in blue and selection in red as in paper
  else
      postcolor = RGBA(0.75, 0.3, 0.3)
  end

  bar(DFreshist[:VAF], DFreshist[:truecounts], linecolor = RGBA(0.431, 0.431, 0.431, 0.9),
  fillcolor = RGBA(0.431, 0.431, 0.431, 0.9))
  plot!(DFres[:VAF], DFres[:mean], color = postcolor, w = 2)
  plot!(DFres[:VAF], DFres[:upperq95], fillrange = DFres[:lowerq95],
               fillalpha = 0.5,
               fillcolor = postcolor, linecolor = false,
               markerstrokecolor=:white, titlefont = font(8, "Calibri"), ytickfont = font(6, "Calibri"), xtickfont = font(6, "Calibri"), legend = false, grid = false,
               yaxis = ("Counts"), xaxis = ("VAF"))
  if model == 2
      xint = median(res.Posterior[model].Parameters[:frequency], Weights(res.Posterior[model].Parameters[:weight])) *
      median(res.Posterior[model].Parameters[:cellularity], Weights(res.Posterior[model].Parameters[:weight]))
      xint = xint/2
      vline!([xint], line=(2,:dash,0.6,:black))
  end

  if model == 3
      xint1 = median(res.Posterior[model].Parameters[:frequency1], Weights(res.Posterior[model].Parameters[:weight])) *
      median(res.Posterior[model].Parameters[:cellularity], Weights(res.Posterior[model].Parameters[:weight]))
      xint1 = xint1 / 2
      xint2 = median(res.Posterior[model].Parameters[:frequency2], Weights(res.Posterior[model].Parameters[:weight])) *
      median(res.Posterior[model].Parameters[:cellularity], Weights(res.Posterior[model].Parameters[:weight]))
      xint12 = xint2 / 2
      vline!([xint1, xint2], line=(2,:dash,0.6,[:black, :black]))
  end

end

function plotmodelposterior(res; plotsbackend = nothing)
    if plotsbackend != nothing
        plotsbackend()
    end
    DF = DataFrame(Model = map(x -> "$x", res.ModelProb[:Model]), Probability = res.ModelProb[:Probability])

    Plots.bar(DF[:Model], DF[:Probability],
    title="Model Probabilities", yaxis = ("Probability"), xaxis = ("Number of subclones"),
    linecolor = :white, fillcolor = RGBA(0.5, 0.5, 0.5, 0.8),
    markerstrokecolor=:white, titlefont = font(10, "Calibri"), ytickfont = font(8, "Calibri"), xtickfont = font(8, "Calibri"), legend = false, grid = false)
end

function plotparameterposterior(res, model = 1; plotsbackend = nothing)
    if plotsbackend != nothing
        plotsbackend()
    end
    if model == 0
        model = model + 1
        Plots.histogram(Array(res.Posterior[model].Parameters)[:, 1:3], nbins = 20, layout = 3, weights = Array(res.Posterior[model].Parameters[:weight]),
        title=["Mutation rate" "# Clonal mutations" "Cellularity"],
        linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
        markerstrokecolor=:white, titlefont = font(10, "Calibri"), ytickfont = font(6, "Calibri"), xtickfont = font(6, "Calibri"), legend = false)
    elseif model == 1
        model = model + 1
        Plots.histogram(Array(res.Posterior[model].Parameters)[:, 1:7],
        nbins = 20, layout = 7,weights = Array(res.Posterior[model].Parameters[:weight]),
        title=["Mutation rate" "# Clonal mutations" "s" "t" "Cellularity" "Subclone frequency" "# Mutations in subclone"],
        linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
        markerstrokecolor=:white, titlefont = font(10, "Calibri"),
        xtickfont = font(6, "Calibri"), ytickfont = font(6, "Calibri"), legend = false)
    elseif model == 2
        model = model + 1
        Plots.histogram(Array(res.Posterior[model].Parameters)[:, 1:11],
        nbins = 20, layout = 11,weights = Array(res.Posterior[model].Parameters[:weight]),
        title=["Mutation rate" "# Clonal mutations" "s1" "t1" "s2" "t2" "Cellularity" "Subclone 1 frequency" "Subclone 2 frequency" "# Mutations in \nsubclone 1" "# Mutations in subclone 2"],
        linecolor = :white, titlefont = font(8, "Calibri"),xtickfont = font(6, "Calibri"), ytickfont = font(6, "Calibri"), fillcolor = RGBA(0.75, 0.3, 0.3), legend = false)
    end
end

"""
    saveallplots(res; <keyword arguments>)

Create and save all plots. For each model, the VAF histogram with model results overlayed is plotted as is the posterior distributions for parameters and models. Function takes a fitABCmodels output type.
...
## Arguments
- `resultsdirectory = "output"`: Directory to save the plots.
- `outputformat = ".pdf"`: Format to save the plots as. Default is pdf, other common formats are available such png etc.
- `plotsbackend = nothing`: Backend to use for plotting in Plots.jl
...
"""
function saveallplots(res; resultsdirectory = "output", outputformat = ".pdf", plotsbackend = nothing)
    if plotsbackend != nothing
        plotsbackend()
    end
  sname = res.SampleName
  dir = joinpath(resultsdirectory, res.SampleName)
  makedirectory(resultsdirectory)
  makeplotsdirectories(dir)
  p = plotmodelposterior(res, plotsbackend = plotsbackend)
  savefig(joinpath(dir, "plots", "$(sname)-modelposterior$(outputformat)"))

  model = 0
  for post in res.Posterior
    if post.Probability > 0.0
      p = plothistogram(res, model, plotsbackend = plotsbackend)
      savefig(joinpath(dir, "plots", "$(sname)-histogram-$(model)clone$(outputformat)"))
      p = plotparameterposterior(res, model, plotsbackend = plotsbackend);
      savefig(joinpath(dir, "plots", "$(sname)-posterior-$(model)clone$(outputformat)"))
    end
    model = model + 1
  end

end
