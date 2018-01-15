using SubClonalSelection

srand(1)
#5509 seconds
@time out = fitABCmodels("example/oneclone.txt",
  "oneclone",
  read_depth = 300,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 2*10^5,
  Nmax = 2^12,
  maxclones = 1,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^6,
  fmin = 0.01);
saveallplots(out, resultsdirectory = "example/")

srand(1)
#5509 seconds
@time out = fitABCmodels("example/oneclone-2.txt",
  "oneclone-2",
  read_depth = 300,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 2*10^5,
  Nmax = 2^12,
  maxclones = 1,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^6,
  fmin = 0.01);
saveallplots(out, resultsdirectory = "example/")


srand(123)
#3379 seconds
@time out = fitABCmodels("example/neutral.txt",
  "neutral",
  read_depth = 300,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 2*10^5,
  maxclones = 1,
  save = true,
  Nmax = 2^12,
  adaptpriors = true,
  Nmaxinf = 10^6,
  fmin = 0.01);
saveallplots(out, resultsdirectory = "example/")


srand(123)
#6068 seconds
@time out = fitABCmodels("example/nikzainal.txt",
  "nikzainal",
  read_depth = 180,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 2*10^5,
  minreads = 6,
  maxclones = 2,
  save = true);
saveallplots(out, resultsdirectory = "example/")


ind=map(x->x.model, out.ABCresults.particles).==3
p = out.ABCresults.particles[ind]
DF = p[10].other[1]

using Gadfly
plot(DF, x = :VAF, y = :freq, Geom.bar,
Theme(default_color = RGBA(0.5, 0.5, 0.5, 0.8)))
