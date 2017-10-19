library(qtl)
library(snowfall)

load("~/F2/output/LG.f2.madmapper.before.crossover.ripple1.Rdata")

sfInit(parallel = TRUE, cpus = 12)
sfLibrary(qtl)
sfExport("LG.f2.madmapper.before.crossover")

rip2 <-
  sfLapply(names(LG.f2.madmapper.before.crossover$geno), function(LG) {
    ripple(cross = LG.f2.madmapper.before.crossover, window = 3, chr = LG, method = "likelihood", error.prob = 0.001)
  })

sfStop()

save(rip2, file = "~/F2/output/missing_rate_0.10/rip2.Rdata")
