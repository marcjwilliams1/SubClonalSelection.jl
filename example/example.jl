using SubClonalSelection

srand(123)
#5509 seconds
@time out = fitABCmodels("example/oneclone.txt",
  "oneclone",
  read_depth = 150,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 3 * 10^5,
  Nmax = 10^3,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true,
  adaptpriors = true,
  Nmaxinf = 10^6);
saveallplots(out, resultsdirectory = "example/")


srand(123)
#3379 seconds
@time out = fitABCmodels("example/neutral.txt",
  "neutral",
  read_depth = 150,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 3 * 10^5,
  Nmax = 10^3,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true,
  adaptpriors = true,
  Nmaxinf = 10^6);
saveallplots(out, resultsdirectory = "example/")


srand(123)
#6068 seconds
@time out = fitABCmodels("example/nikzainal.txt",
  "nikzainal",
  read_depth = 180,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 3 * 10^5,
  minreads = 6,
  Nmax = 10^3,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true);
saveallplots(out, resultsdirectory = "example/")


ind=map(x->x.model, out.ABCresults.particles).==3
p = out.ABCresults.particles[ind]
DF = p[10].other[1]

using Gadfly
plot(DF, x = :VAF, y = :freq, Geom.bar,
Theme(default_color = RGBA(0.5, 0.5, 0.5, 0.8)))
