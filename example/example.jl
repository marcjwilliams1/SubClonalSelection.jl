using SubClonalSelection

srand(123)
#5509 seconds
@time out = fitABCmodels("example/oneclone.txt",
  "oneclone",
  read_depth = 300,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 10^4,
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
  maxiterations = 10^4,
  maxclones = 1,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^6,
  fmin = 0.01);
saveallplots(out, resultsdirectory = "example/")

srand(123)
@time out = fitABCmodels("example/4990-12/data/4990-12.txt",
  "4990-12",
  read_depth = 150,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 10^4,
  maxclones = 1,
  ρ = 0.005,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^10,
  minvaf = 0.04,
  mincellularity = 0.95,
  fmin = 0.04);
saveallplots(out, resultsdirectory = "example/otherdata/")


srand(1)
#5509 seconds
@time out = fitABCmodels("example/oneclone-2.txt",
  "oneclone-2",
  read_depth = 300,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 5*10^5,
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
@time out = fitABCmodels("example/oneclone-3.txt",
  "oneclone-3",
  read_depth = 500,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 5*10^5,
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
@time out = fitABCmodels("example/oneclone-4.txt",
  "oneclone-4",
  read_depth = 300,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 5*10^5,
  Nmax = 2^12,
  maxclones = 1,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^0,
  fmin = 0.01);
saveallplots(out, resultsdirectory = "example/")



srand(1)
@time out = fitABCmodels("example/otherdata/data/MO_1444.txt",
  "MO_1044",
  read_depth = 171,
  resultsdirectory = "example/otherdata/",
  nparticles = 100,
  maxiterations = 10^5,
  Nmax = 2^12,
  maxclones = 1,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^10,
  minreads = round(Int64, 171*0.05),
  fmin = 0.05);
saveallplots(out, resultsdirectory = "example/otherdata/")

srand(1)
@time out = fitABCmodels("example/otherdata/data/MO_1031.txt",
  "MO_1031",
  read_depth = 151,
  resultsdirectory = "example/otherdata/",
  nparticles = 100,
  maxiterations = 10^5,
  Nmax = 2^12,
  maxclones = 1,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^10,
  minreads = round(Int64, 151*0.05),
  fmin = 0.05);
saveallplots(out, resultsdirectory = "example/otherdata/")

srand(1)
@time out = fitABCmodels("example/otherdata/data/MO_1425.txt",
  "MO_1425",
  read_depth = 185,
  resultsdirectory = "example/otherdata/",
  nparticles = 100,
  maxiterations = 10^5,
  Nmax = 2^12,
  maxclones = 1,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^10,
  minreads = round(Int64, 185*0.05),
  fmin = 0.05);
saveallplots(out, resultsdirectory = "example/otherdata/")

srand(123)
@time out = fitABCmodels("example/otherdata/data/4990-12.txt",
  "4990-12",
  read_depth = 150,
  resultsdirectory = "example/otherdata/",
  nparticles = 100,
  maxiterations = 10^5,
  maxclones = 2,
  ρ = 0.005,
  save = true,
  adaptpriors = true,
  verbose = true,
  Nmaxinf = 10^10,
  minvaf = 0.04,
  mincellularity = 0.95,
  fmin = 0.04);
saveallplots(out, resultsdirectory = "example/otherdata/")

srand(1)
@time out = fitABCmodels("example/otherdata/data/4990-14.txt",
  "4990-14",
  read_depth = 150,
  resultsdirectory = "example/otherdata/",
  nparticles = 100,
  maxiterations = 10^5,
  Nmax = 2^12,
  maxclones = 1,
  ρ = 0.005,
  save = true,
  adaptpriors = false,
  verbose = true,
  Nmaxinf = 10^10,
  minreads = round(Int64, 150*0.05),
  mincellularity = 0.95,
  fmin = 0.05);
saveallplots(out, resultsdirectory = "example/otherdata/")


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
  verbose = true,
  Nmaxinf = 10^10,
  fmin = 0.05);
saveallplots(out, resultsdirectory = "example/")


srand(123)
#6068 seconds
@time out = fitABCmodels("example/nikzainal.txt",
  "nikzainal",
  read_depth = 180,
  resultsdirectory = "example/",
  nparticles = 100,
  maxiterations = 2*10^5,
  minvaf = 0.04,
  adaptpriors = true,
  maxclones = 2,
  save = true);
saveallplots(out, resultsdirectory = "example/")
