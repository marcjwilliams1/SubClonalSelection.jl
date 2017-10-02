using SubClonalSelection

@time out = fitABCmodels("/Users/marcwilliams/Google\ Drive/early_events/Revisions/methodvalidation/2clonewgs/1CloneWGS.larget-7.txt",
  "1CloneWGS.larget-7",
  read_depth = 150,
  resultsdirectory = "/Users/marcwilliams/Google\ Drive/early_events/Revisions/methodvalidation/tests/",
  nparticles = 250,
  maxiterations = 7*10^4,
  Nmax = 10^3,
  Nmaxinf = 10^6,
  maxclones = 2,
  save = true,
  firstpass = false,
  verbose = true,
  adaptpriors = true);

saveallplots(out, resultsdirectory = "/Users/marcwilliams/Google\ Drive/early_events/Revisions/methodvalidation/tests/")
