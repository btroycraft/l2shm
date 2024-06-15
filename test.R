install.packages("./", repos=NULL, type="source")
library(l2shm)

U <- {
  data <- read.csv("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv", header=TRUE)
  data <- subset(data, location == 'Westchester Lagoon', c(-location, -sample))
  u.scores(data)
}

par <- l2shm.gd.emp(U, .3^2, terms = 35L, maxiter = 10, groups = 100L, tol=-1)
par <- l2shm.gd.emp(U, .15^2, terms = 35L, maxiter = 200, par = par, tol=-1)

par_sub_list <- lapply(1:groups, function(groups_sub){
  l2shm.gd.proj(par, groups = groups_sub, maxiter = 10^3)
})
nrm2 <- sapply(par_sub_list, function(par_sub){
  l2shm.nrm2.sq.diff(par_sub, par, terms=50L)
})
nrm2 <- sapply(1:groups, function(ind) min(nrm2[1:ind]))
nrm2[groups] <- 0

library(ggplot2)

t <- seq(0, 2, .001)
l <- sapply(t, function(t){
  min(which(nrm2<.41+t))
})
u <- sapply(t, function(t){
  max(which(nrm2+.41>t))
})

plot(t, u, type='l')
lines(t, l)
qplot(t, u, t, l)
qplot(t, l, add=TRUE)

ggplot() +
  geom_line(aes(x = t, y = l), col='red') +
  geom_line(aes(x = t, y = u), col='blue') +
  xlab(expression(epsilon)) +
  ylab(expression(kappa(epsilon)))

t <- seq(-1, 1, .01); plot(t, heat.kern(t, .1^2, nrow(U), 35L), type='l')

{
  set.seed(5675867)

  groups0 <- 4
  t0 <- .1^2
  par0 <- list(
    mu = runif.sph(groups0, 3),
    t = rep(t0, groups0),
    alpha = local((._ <- runif(groups0))/sum(._))
  )

  groups <- 50
  tmin <- .05^2

  library(parallel)
  RNGkind("L'Ecuyer-CMRG")

  set.seed(5417815)

  CORES <- 6

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
    unlist(mclapply(1:(ceiling(10^2/CORES)*CORES), function(._){
      runif(1)
      U_boot <- rheat.sph.mix(ncol(U), par_sub)
      par_boot <- l2shm.gd.emp(U_boot, tmin, maxiter = 10^3, groups = groups, tol = 10^-7)
      l2shm.nrm2.sq.diff(par_boot, par_sub)
    }, mc.cores = CORES))
  }, simplify = FALSE)


  plot.new()
  plot(c(), xlim=c(0, 1), ylim=c(0, 1), xlab='Nominal', ylab='Actual')
  segments(0, 0, 1, 1)

  t <- seq(0, 1, .001)
  for(reps_boot in reps_boot_list){
    lines(t, cdf(reps)(icdf(reps_boot)(t)), xlim=c(0, 1), ylim=c(0, 1))
  }

  plot(t, cdf(reps)(t), ylim=c(0, 1), xlim=c(0, 1), type='l'); lines(t, cdf(reps_boot)(t))

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

  t <- seq(0, 1, length.out=1000)
  plot(t, cdf(reps)(icdf(reps_boot)(t)), xlim=c(0, 1), ylim=c(0, 1), type='l', asp=1)
  abline(0, 1)

  library(plotly)
  U <- rheat.sph.mix(10^3, par0)
  plot_ly() %>%
    add_trace(x = U[1, ], y = U[2, ], z = U[3, ], mode = "markers", marker=list(size=2)) %>%
    add_trace(x = 1.1*par$mu[1, ], y = 1.1*par$mu[2, ], z = 1.1*par$mu[3, ], mode = "markers", color=par$alpha) #%>%
    #add_trace(x = 1.2*par_sub$mu[1, ], y = 1.2*par_sub$mu[2, ], z = 1.2*par_sub$mu[3, ], mode = "markers", color=par_sub$alpha)

  U_boot <- rheat.sph.mix(10^3, par)

  plot_ly() %>%
    add_trace(x = U[1, ], y = U[2, ], z = U[3, ], mode = "markers", marker=list(size=2)) %>%
    add_trace(x = U_boot[1, ], y = U_boot[2, ], z = U_boot[3, ], mode = "markers", marker=list(size=2)) %>%
    add_trace(x = 1.1*par$mu[1, ], y = 1.1*par$mu[2, ], z = 1.1*par$mu[3, ], mode = "markers", color=par$alpha)

  l2.heat.obj.emp(U, y)
}
