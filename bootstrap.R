install.packages("./", repos=NULL, type="source")
library(l2shm)
library(parallel)

cdf <- function(x){
  x <- sort(x)
  function(y) sapply(y, function(y){
    if(y <= x[1]){
      0
    } else if(y >= x[length(x)]) {
      1
    } else {
      i_low <- max(which(y >= x))
      i_upp <- min(which(y <= x))
      (i_low + (y-x[i_low])/(x[i_upp]-x[i_low]) - 1) / (length(x)-1)
    }
  })
}

icdf <- function(x){
  x <- sort(x)
  function(y) sapply(y, function(y){
    if(y == 0){
      x[1]
    } else if(y == 1){
      x[length(x)]
    } else {
      i_low = floor(y*(length(x)-1))
      x[i_low+1] + (y*(length(x)-1)-i_low)*(x[i_low+2]-x[i_low+1])
    }
  })
}

RNGkind("L'Ecuyer-CMRG")
set.seed(5675867)

CORES <- 6

groups0 <- 4
t0 <- .1^2
par0 <- list(
  mu = runif.sph(groups0, 3),
  t = rep(t0, groups0),
  alpha = local((._ <- runif(groups0))/sum(._))
)

groups <- 50
tmin <- .05^2


reps <- unlist(mclapply(1:(ceiling(10^2/CORES)*CORES), function(ind){
  runif(1)
  U <- rheat.sph.mix(10^3, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)
  l2shm.nrm2.sq.diff(par, par0)
}, mc.cores = CORES))

reps_boot_list <- replicate(10, {
  U <- rheat.sph.mix(10^3, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^3, groups = groups, tol = 10^-7)
  par_sub <- {
    par_sub_list <- lapply(1:groups, function(groups_sub){
      l2shm.gd.proj(par, groups = groups_sub, maxiter = 10^3)
    })
    par_sub_list[[groups]] <- copy.par(par)
    nrm2 <- sapply(par_sub_list, function(par_sub){
      l2shm.nrm2.sq.diff(par_sub, par, terms=50L)
    })
    nrm2[groups] <- 0;
    par_sub_list[[min(which(nrm2 <= 44/ncol(U)^.7))]]
  }
  unlist(mclapply(1:(ceiling(3*10^2/CORES)*CORES), function(._){
    runif(1)
    U_boot <- rheat.sph.mix(ncol(U), par_sub)
    par_boot <- l2shm.gd.emp(U_boot, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)
    l2shm.nrm2.sq.diff(par_boot, par_sub)
  }, mc.cores = CORES))
}, simplify = FALSE)

{
  plot.new()
  plot(c(),
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab='Nominal',
       ylab='Actual'
  )
  segments(0, 0, 1, 1)

  t <- seq(0, 1, .001)
  for(reps_boot in reps_boot_list){
    lines(t, cdf(reps)(icdf(reps_boot)(t)), xlim=c(0, 1), ylim=c(0, 1))
  }
}

save(reps, reps_boot_list, file = "save.RData")

load("save.RData")
