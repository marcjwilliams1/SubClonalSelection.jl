function plothistogram(res, model = 0)

  model = model + 1
  DFres = res.Posterior[model].MeanHistogram
  DFres[:VAF] = collect(0.01:0.01:1.0)
  DF = DataFrame(VAF = res.VAF)
  DFres = DFres[1:75, :]

  l1 = layer(DFres, x = :VAF, y = :mean, ymin = :lowerq95, ymax = :upperq95, Geom.line, Geom.ribbon,
  Theme(line_width = 0.06cm, default_color = RGBA(0.75, 0.3, 0.3),
  lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.5)))
  l2 = layer(DFres, x = :VAF, y = :truecounts, Geom.bar,
  Theme(default_color = RGBA(0.5, 0.5, 0.5, 0.8),
  major_label_font_size = 16pt,
  minor_label_font_size = 12pt))

  myplot = Gadfly.plot(l1, l2,
  Guide.xlabel("VAF"),
  Guide.ylabel("Counts"))

  return myplot
end

function plotmodelposterior(res)
    p = Gadfly.plot(res.ModelProb, x=:Model, y = :Probability, Geom.bar,
    Theme(bar_spacing = 0.2cm,
    default_color = RGBA(0.5, 0.5, 0.5, 0.8),
    major_label_font_size = 16pt,
    minor_label_font_size = 12pt))

    return p
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
  draw(PNG(joinpath(dir, "plots", "$(sname)-modelposterior.png"), 4inch, 3inch), p)

  model = 0
  for post in res.Posterior
    if post.Probability > 0.0
      p = plothistogram(res, model)
      draw(PNG(joinpath(dir, "plots", "$(sname)-histogram-$(model)clone.png"), 4inch, 3inch), p)
      p = plotparameterposterior(res, model);
      savefig(p, joinpath(dir, "plots", "$(sname)-posterior-$(model)clone.png"))
    end
    model = model + 1
  end

end
