function plothistogram(res, model = 0)

  model = model + 1
  DFres = res.Posterior[model].MeanHistogram
  DFres[:VAF] = collect(0.01:0.01:1.0)
  DF = DataFrame(VAF = res.VAF)

  l1 = layer(DFres, x = :VAF, y = :mean, ymin = :lowerq, ymax = :upperq, Geom.line, Geom.ribbon,
  Theme(default_color = RGBA(0.75, 0.3, 0.3),
  lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.5)))
  l2 = layer(DF, x = :VAF, Geom.histogram(bincount=100),
  Theme(default_color = RGBA(0.5, 0.5, 0.5, 0.8),
  major_label_font_size = 16pt,
  minor_label_font_size = 12pt))

  myplot = plot(l1, l2,
  Guide.xlabel("VAF"),
  Guide.ylabel("Counts"))

  return myplot
end

function plotmodelposterior(res)
    p = plot(res.ModelProb, x=:Model, y = :Probability, Geom.bar,
    Theme(bar_spacing = 0.2cm,
    default_color = RGBA(0.5, 0.5, 0.5, 0.8),
    major_label_font_size = 16pt,
    minor_label_font_size = 12pt))

    return p
end

function plotparameterposterior(res, model = 1)

    model = model + 1

    DF = stack(res.Posterior[model].Parameters)
    plot(DF, x=:value, xgroup=:variable,
    Geom.subplot_grid(Geom.histogram(bincount = 50), free_x_axis=true),
    Theme(
    default_color = RGBA(0.5, 0.5, 0.5, 0.8),
    major_label_font_size = 12pt,
    minor_label_font_size = 8pt))

end

function saveallplots(res; resultsdirectory = "output")

  makeplotsdirectories(resultsdirectory)
  p = plotmodelposterior(res)
  draw(PNG(joinpath(resultsdirectory, "plots", "modelposterior.png"), 4inch, 3inch), p)

  model = 0
  for post in res.Posterior
    if post.Probability > 0.0
      p = plothistogram(res, model + 1)
      draw(PNG(joinpath(resultsdirectory, "plots", "histogram-$(model)clone.png"), 4inch, 3inch), p)
      p = plotparameterposterior(res, model + 1)
      draw(PNG(joinpath(resultsdirectory, "plots", "posterior-$(model)clone.png"), 15inch, 6inch), p)
    end
    model = model + 1
  end

end
