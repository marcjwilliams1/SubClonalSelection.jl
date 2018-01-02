using SubClonalSelection

@time out = fitABCmodels("example/oneclone.txt",
  "oneclone",
  read_depth = 150,
  resultsdirectory = "example/",
  nparticles = 200,
  maxiterations = 10^4,
  Nmax = 10^3,
  maxclones = 1,
  save = true,
  firstpass = false,
  verbose = true,
  adaptpriors = true);
saveallplots(out, resultsdirectory = "example/")
