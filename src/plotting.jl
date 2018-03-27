function plothistogram(res, model = 0; removelowfrequencies = true)

  if removelowfrequencies == true
      dl = res.ABCsetup.Models[1].constants[8]
  else
      dl = 0.0
  end

  model = model + 1

  DFres = res.Posterior[model].MeanHistogram
  DFres[:VAF] = collect(0.01:0.01:1.0)
  DF = DataFrame(VAF = res.VAF)
  DFres = DFres[1:75, :]
  DFres = DFres[DFres[:VAF].>=dl, :]

  if model - 1 == 0
      postcolor = RGBA(0.0, 76/255, 153/255) # plot neutral in blue and selection in red as in paper
  else
      postcolor = RGBA(0.75, 0.3, 0.3)
  end

  bar(DFres[:VAF], DFres[:truecounts], linecolor = RGBA(0.431, 0.431, 0.431, 0.9),
  fillcolor = RGBA(0.431, 0.431, 0.431, 0.9))
  plot!(DFres[:VAF], DFres[:mean], color = postcolor, w = 2)
  plot!(DFres[:VAF], DFres[:upperq95], fillrange = DFres[:lowerq95],
               fillalpha = 0.5,
               fillcolor = postcolor, linecolor = false,
               markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false, grid = false,
               yaxis = ("Counts"), xaxis = ("VAF"))

end

function plotmodelposterior(res)

    DF = DataFrame(Model = map(x -> "$x", res.ModelProb[:Model]), Probability = res.ModelProb[:Probability])

    Plots.bar(DF[:Model], DF[:Probability],
    title="Model Probabilities", yaxis = ("Probability"), xaxis = ("# Number of subclones"),
    linecolor = :white, fillcolor = RGBA(0.5, 0.5, 0.5, 0.8),
    markerstrokecolor=:white, titlefont = font(14, "Calibri"), ytickfont = font(12, "Calibri"), xtickfont = font(12, "Calibri"), legend = false, grid = false)
end

function plotparameterposterior(res, model = 1)

    if model == 0
        model = model + 1
        Plots.histogram(Array(res.Posterior[model].Parameters)[:, 1:3], nbins = 20, layout = 3,
        title=["Mutation rate" "# Clonal mutations" "Cellularity"],
        linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
        markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false)
    elseif model == 1
        model = model + 1
        Plots.histogram(Array(res.Posterior[model].Parameters)[:, 1:7],
        nbins = 20, layout = 7,
        title=["Mutation rate" "# Clonal mutations" "s" "t" "Cellularity" "Subclone frequency" "# Mutations in subclone"],
        linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
        markerstrokecolor=:white, titlefont = font(12, "Calibri"),
        xtickfont = font(10, "Calibri"), ytickfont = font(10, "Calibri"), legend = false)
    elseif model == 2
        model = model + 1
        Plots.histogram(Array(res.Posterior[model].Parameters)[:, 1:11],
        nbins = 20, layout = 11,
        title=["Mutation rate" "# Clonal mutations" "s1" "t1" "s2" "t2" "Cellularity" "Subclone 1 frequency" "Subclone 2 frequency" "# Mutations in \nsubclone 1" "# Mutations in subclone 2"],
        linecolor = :white, titlefont = font(12, "Calibri"),xtickfont = font(10, "Calibri"), ytickfont = font(10, "Calibri"), fillcolor = RGBA(0.75, 0.3, 0.3), legend = false)
    end

end

function saveallplots(res; resultsdirectory = "output")

  sname = res.SampleName
  dir = joinpath(resultsdirectory, res.SampleName)
  makedirectory(resultsdirectory)
  makeplotsdirectories(dir)
  p = plotmodelposterior(res)
  savefig(joinpath(dir, "plots", "$(sname)-modelposterior.png"))

  model = 0
  for post in res.Posterior
    if post.Probability > 0.0
      p = plothistogram(res, model)
      savefig(joinpath(dir, "plots", "$(sname)-histogram-$(model)clone.png"))
      p = plotparameterposterior(res, model);
      savefig(joinpath(dir, "plots", "$(sname)-posterior-$(model)clone.png"))
    end
    model = model + 1
  end

end
