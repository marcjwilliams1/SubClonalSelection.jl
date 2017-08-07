using SubClonalSelection

out = fitABCmodels("example/neutral.txt",
  "neutral",
  read_depth = 150,
  resultsdirectory = "/Users/marcwilliams/Google\ Drive/test/",
  nparticles = 200,
  maxiterations = 5*10^3,
  Nmax = 10^3,
  maxclones = 1,
  save = true,
  firstpass = false);
