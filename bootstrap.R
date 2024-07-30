install.packages("./", repos=NULL, type="source")
library(l2shm)
library(parallel)
library(ggplot2)

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

CORES <- 4

setDefaultCluster(makeCluster(CORES))
invisible(clusterEvalQ(getDefaultCluster(), library(l2shm)))

{
  groups0 <- 4
  t0 <- .1^2
  par0 <- list(
    mu = runif.sph(groups0, 3),
    t = rep(t0, groups0),
    alpha = local((._ <- runif(groups0))/sum(._))
  )

  groups <- 50
  tmin <- .05^2
}

{
  par_list0 <- local({
    par_list <- lapply(1:(groups0-1), function(groups_proj){
      l2shm.gd.proj(par0, groups = groups_proj, maxiter = 10^5, terms = 100L)
    })
    for(ind in groups0:groups){
      par_list[[ind]] <- copy.par(par0)
    }
    par_list
  })
  nrm0 <- local({
    nrm <- sapply(par_list0[1:(groups0-1)], function(par_proj){
      sqrt(l2shm.nrm2.sq.diff(par_proj, par0, terms = 100L))
    })
    nrm <- c(sqrt(l2shm.nrm2.sq(par0, terms = 100L)-1), nrm, rep(0, groups-groups0+1))
    nrm <- sapply(seq_along(nrm), function(ind) min(nrm[1:ind]))
    cbind(k = seq.int(0, groups), n = nrm)
  })
}

system.time({
  U <- rheat.sph.mix(10^4, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^3, groups = groups, tol = 10^-7, terms = 100L)
  par_list <- lapply(1:(groups-1), function(groups_proj){
    l2shm.gd.proj(par, groups = groups_proj, maxiter = 10^3, terms = 100L)
  })
  nrm <- local({
    nrm1 <- sapply(par_list[1:(groups-1)], function(par_proj){
      sqrt(l2shm.nrm2.sq.diff(par_proj, par, terms = 100L))
    })
    nrm1 <- c(sqrt(l2shm.nrm2.sq(par, terms = 100L)-1), nrm1, 0)
    nrm1 <- sapply(seq_along(nrm1), function(ind) min(nrm1[1:ind]))
    nrm2 <- sapply(par_list0[1:groups0], function(par_proj){
      sqrt(l2shm.nrm2.sq.diff(par_proj, par, terms = 100L))
    })
    nrm2 <- c(sqrt(l2shm.nrm2.sq(par, terms = 100L)-1), nrm2, rep(nrm2[groups0], groups-groups0))
    nrm2 <- sapply(seq_along(nrm2), function(ind) min(nrm2[1:ind]))
    cbind(k = seq.int(0, groups), n = pmin(nrm1, nrm2))
  })
})

clusterExport(getDefaultCluster(), c("nrm0", "groups", "tmin", "groups0", "t0", "par0", "par_list0"), envir = environment())


reps_cov_k4_p1000 <- unlist(parLapply(1:(ceiling(10^3/CORES)*CORES), function(ind){
  runif(1)
  U <- rheat.sph.mix(10^3, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^3, groups = groups, tol = 10^-7, terms = 100L)
  par_list <- lapply(1:(groups-1), function(groups_proj){
    l2shm.gd.proj(par, groups = groups_proj, maxiter = 10^3, terms = 100L)
  })
  nrm <- local({
    nrm1 <- sapply(par_list[1:(groups-1)], function(par_proj){
      sqrt(l2shm.nrm2.sq.diff(par_proj, par, terms = 100L))
    })
    nrm1 <- c(sqrt(l2shm.nrm2.sq(par, terms = 100L)-1), nrm1, 0)
    nrm1 <- sapply(seq_along(nrm1), function(ind) min(nrm1[1:ind]))
    nrm2 <- sapply(par_list0[1:groups0], function(par_proj){
      sqrt(l2shm.nrm2.sq.diff(par_proj, par, terms = 100L))
    })
    nrm2 <- c(sqrt(l2shm.nrm2.sq(par, terms = 100L)-1), nrm2, rep(nrm2[groups0], groups-groups0))
    nrm2 <- sapply(seq_along(nrm2), function(ind) min(nrm2[1:ind]))
    cbind(k = seq.int(0, groups), n = pmin(nrm1, nrm2))
  })
  max(abs(nrm[, "n"]-nrm0[, "n"]))
}, cl = getDefaultCluster()))

save(reps_cov_k4_p1000, file = "reps_cov_k4_p1000.RData")

reps_cov_k4_p10000 <- unlist(parLapply(1:(ceiling(10^3/CORES)*CORES), function(ind){
  runif(1)
  U <- rheat.sph.mix(10^4, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^3, groups = groups, tol = 10^-7, terms = 100L)
  par_list <- lapply(1:(groups-1), function(groups_proj){
    l2shm.gd.proj(par, groups = groups_proj, maxiter = 10^3, terms = 100L)
  })
  nrm <- local({
    nrm1 <- sapply(par_list[1:(groups-1)], function(par_proj){
      sqrt(l2shm.nrm2.sq.diff(par_proj, par, terms = 100L))
    })
    nrm1 <- c(sqrt(l2shm.nrm2.sq(par, terms = 100L)-1), nrm1, 0)
    nrm1 <- sapply(seq_along(nrm1), function(ind) min(nrm1[1:ind]))
    nrm2 <- sapply(par_list0[1:groups0], function(par_proj){
      sqrt(l2shm.nrm2.sq.diff(par_proj, par, terms = 100L))
    })
    nrm2 <- c(sqrt(l2shm.nrm2.sq(par, terms = 100L)-1), nrm2, rep(nrm2[groups0], groups-groups0))
    nrm2 <- sapply(seq_along(nrm2), function(ind) min(nrm2[1:ind]))
    cbind(k = seq.int(0, groups), n = pmin(nrm1, nrm2))
  })
  max(abs(nrm[, "n"]-nrm0[, "n"]))
}, cl = getDefaultCluster()))

save(reps_cov_k4_p10000, file = "reps_cov_k4_p10000.RData")


reps_nrm_k4_p1000 <- unlist(parLapply(1:(ceiling(10^3/CORES)*CORES), function(ind){
  runif(1)
  U <- rheat.sph.mix(10^3, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)

  sqrt(l2shm.nrm2.sq.diff(par, par0))
}, cl = getDefaultCluster()))

save(reps_nrm_k4_p1000, file = "reps_nrm_k4_p1000.RData")

reps_nrm_k4_p10000 <- unlist(parLapply(1:(ceiling(10^2/CORES)*CORES), function(ind){
  runif(1)
  U <- rheat.sph.mix(10^4, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)

  sqrt(l2shm.nrm2.sq.diff(par, par0))
}, cl = getDefaultCluster()))

save(reps_nrm_k4_p10000, file = "reps_nrm_k4_p10000.RData")


reps_boot_list_k4_p1000 <- replicate(10, {
  U <- rheat.sph.mix(10^3, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^3, groups = groups, tol = 10^-7)
  clusterExport(getDefaultCluster(), c("U", "par"), envir = environment())
  unlist(parLapply(1:(ceiling(3*10^2/CORES)*CORES), function(._){
    runif(1)
    U_boot <- rheat.sph.mix(ncol(U), par)
    par_boot <- l2shm.gd.emp(U_boot, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)
    sqrt(l2shm.nrm2.sq.diff(par_boot, par))
  }, cl = getDefaultCluster()))
}, simplify = FALSE)

save(reps_boot_list_k4_p1000, file = "reps_boot_list_k4_p1000.RData")

reps_boot_list_k4_p10000 <- replicate(10, {
  U <- rheat.sph.mix(10^4, par0)
  par <- l2shm.gd.emp(U, tmin, maxiter = 10^3, groups = groups, tol = 10^-7)
  clusterExport(getDefaultCluster(), c("U", "par"), envir = environment())
  unlist(parLapply(1:(ceiling(3*10^2/CORES)*CORES), function(._){
    runif(1)
    U_boot <- rheat.sph.mix(ncol(U), par)
    par_boot <- l2shm.gd.emp(U_boot, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)
    sqrt(l2shm.nrm2.sq.diff(par_boot, par))
  }, cl = getDefaultCluster()))
}, simplify = FALSE)

save(reps_boot_list_k4_p10000, file = "reps_boot_list_k4_p10000.RData")


reps_boot_list_k4_p1000_trans <- replicate(10, {
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
    par_sub_list[[sum(nrm2 > 44/ncol(U)^.7)+1]]
  }
  clusterExport(getDefaultCluster(), c("U", "par_sub"), envir = environment())
  unlist(parLapply(1:(ceiling(3*10^2/CORES)*CORES), function(._){
    runif(1)
    U_boot <- rheat.sph.mix(ncol(U), par_sub)
    par_boot <- l2shm.gd.emp(U_boot, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)
    sqrt(l2shm.nrm2.sq.diff(par_boot, par_sub))
  }, cl = getDefaultCluster()))
}, simplify = FALSE)

save(reps_boot_list_k4_p1000_trans, file = "reps_boot_list_k4_p1000_trans.RData")

reps_boot_list_k4_p10000_trans <- replicate(10, {
  U <- rheat.sph.mix(10^4, par0)
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
    par_sub_list[[sum(nrm2 > 44/ncol(U)^.7)+1]]
  }
  clusterExport(getDefaultCluster(), c("U", "par_sub"), envir = environment())
  unlist(parLapply(1:(ceiling(3*10^2/CORES)*CORES), function(._){
    runif(1)
    U_boot <- rheat.sph.mix(ncol(U), par_sub)
    par_boot <- l2shm.gd.emp(U_boot, tmin, maxiter = 10^2, groups = groups, tol = 10^-7)
    sqrt(l2shm.nrm2.sq.diff(par_boot, par_sub))
  }, cl = getDefaultCluster()))
}, simplify = FALSE)

save(reps_boot_list_k4_p10000_trans, file = "reps_boot_list_k4_p10000_trans.RData")



load("reps_nrm_k4_p1000.RData")
load("reps_nrm_k4_p10000.RData")
load("reps_cov_k4_p1000.RData")
load("reps_cov_k4_p10000.RData")

{
  library(cobs)
  t <- seq(0, 1, .001)

  data <- rbind(
    data.frame(
      nominal = t,
      actual = cdf(reps_nrm_k4_p1000)(icdf(reps_cov_k4_p1000)(t)),
      p = 1000
    ),
    data.frame(
      nominal = t,
      actual = cdf(reps_nrm_k4_p10000)(icdf(reps_cov_k4_p10000)(t)),
      p = 10000
    )
  )
  data$p <-  as.factor(data$p)
  data$actual[data$nominal == 0] <- 0
  data$actual[data$nominal == 1] <- 1
  data$actual <- pmax(data$nominal, data$actual)


  # is <- interpSplineCon(conreg(data$nominal[data$p == 1000], data$actual[data$p == 1000], maxit=c(1000, 100)))
  # data$actual[data$p == 1000] <- predict(is, data$nominal[data$p == 1000])$y
  # data$actual[data$nominal == 0] <- 0
  # data$actual[data$nominal == 1] <- 1
  #
  # is <- interpSplineCon(conreg(data$nominal[data$p == 10000], data$actual[data$p == 10000], maxit=c(1000, 100)))
  # data$actual[data$p == 10000] <- predict(is, data$nominal[data$p == 10000])$y
  # data$actual[data$nominal == 0] <- 0
  # data$actual[data$nominal == 1] <- 1

  plot2 <- ggplot(data = data) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, linewidth = 1) +
    geom_line(aes(x = nominal, y = actual, col = p), linewidth = .8)+
    scale_color_hue(l=50)

  plot2_zoom <- ggplot(data = data) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, linewidth = 1) +
    geom_line(aes(x = nominal, y = actual, col = p), linewidth = .8)+
    scale_color_hue(l=50) +
    coord_cartesian(xlim=c(.9, 1))

  ggsave("quantiles1.png", plot2, scale = 1.3, width = 3, height = 2, units = "in")
}

load("reps_boot_list_k4_p1000.RData")
load("reps_boot_list_k4_p10000.RData")

{
  t <- seq(0, 1, .01)

  data <- rbind(
    do.call(rbind, lapply(seq_along(reps_boot_list_k4_p1000), function(ind){
      data.frame(
        nominal = t,
        actual = cdf(reps_nrm_k4_p1000)(icdf(reps_boot_list_k4_p1000[[ind]])(t)),
        p = 1000,
        run = ind
      )
    })),
    do.call(rbind, lapply(seq_along(reps_boot_list_k4_p10000), function(ind){
      data.frame(
        nominal = t,
        actual = cdf(reps_nrm_k4_p10000)(icdf(reps_boot_list_k4_p10000[[ind]])(t)),
        p = 10000,
        run = ind
      )
    }))
  )

  plot_boot1 <- ggplot(data = data) +
    #geom_line(data = subset(data, p == 1000), aes(x = nominal, y = actual, group = run), col = "darkblue", linewidth = .4) +
    geom_line(data = subset(data, p == 10000), aes(x = nominal, y = actual, group = run), col = "darkred", linewidth = .4) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1)

  ggsave("bootstrap1.png", plot_boot1, scale = 1, width = 3, height = 2, units = "in")
}

load("reps_boot_list2.RData")

{
  t <- seq(0, 1, .01)

  data <- do.call(rbind, lapply(seq_along(reps_boot_list2), function(ind){
    data.frame(
      nominal = t,
      actual = cdf(reps)(icdf(reps_boot_list2[[ind]])(t)),
      run = ind
    )
  }))

  plot_boot2 <- ggplot(data = data) +
    geom_line(aes(x = nominal, y = actual, group = run), col = "darkred", linewidth = .4) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1)

  ggsave("bootstrap2.png", plot_boot2, scale = 1, width = 3, height = 2, units = "in")
}
