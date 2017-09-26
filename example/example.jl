using SubClonalSelection

@time out = fitABCmodels("example/oneclone.txt",
  "oneclone",
  read_depth = 150,
  resultsdirectory = "/Users/marcwilliams/Google\ Drive/test/",
  nparticles = 500,
  maxiterations = 5*10^4,
  Nmax = 10^3,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true,
  adaptpriors = true);

saveallplots(out, resultsdirectory = "/Users/marcwilliams/Google\ Drive/test/")
