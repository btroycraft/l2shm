install.packages("./", repos=NULL, type="source")
library(l2shm)

U <- {
  data <- read.csv("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv", header=TRUE)
  data <- subset(data, location == 'Westchester Lagoon', c(-location, -sample))
  u.scores(data)
}

groups <- 300
tmin <- .05
par <- list(
  mu = unname(U[, sample(ncol(U), groups, replace=TRUE)]),
  t = rep(tmin, groups),
  alpha = rep(1/groups, groups)
)
l2shm.gd.emp(U, tmin, maxiter = 30, par=par, terms=90)
l2.heat(U, par)

l2
y <- heatkern.emp.gd(U, tmin = tmin, start = y, maxiter = 40)

y_max <- heatkern.emp.sgd(U, .02, 4, 10, maxiter = 10)

for(._ in 1:3){
  y <- heatkern.emp.sgd(U, tmin, maxiter = 30, start=y, batch = 100)
  print(l2.heat(U, y)/2)
  y <- heatkern.emp.gd(U, tmin, maxiter = 4, start=y)
  print(l2.heat(U, y)/2)
}




plot(seq(.1, 1, .01)[-1], heat.kern.dt(.9, seq(.1, 1, .01)[-1], 23, 200))
lines(seq(.1, 1, .01)[-1], diff(heat.kern(.9, seq(.1, 1, .01), 23, 200))/diff(seq(.1, 1, .01)))


t <- head(seq(0, 1, .01), -1)

heat.kern.test(.1, 10, 100)
.Call("_heat_kern_test", PACKAGE="l2shm", .1, 10, 100)

{
  groups0 <- 4
  par0 <- list(
    mu = replicate(groups0, (._ <- rnorm(3))/sqrt(sum(._^2))),
    t = rep(.01, groups0),
    alpha = local((._ <- runif(groups0))/sum(._))
  )

  U <- rheat.sph.mix(10^3, par0)

  groups <- 30
  tmin <- .005
  par <- list(
    mu = U[, sample(ncol(U), groups, replace=TRUE)],
    t = rep(2*tmin, groups),
    alpha = rep(1/groups, groups)
  )
  par_cpy <- par

  par <- l2shm.gd.emp(U, tmin, maxiter = 100, groups = groups)
  l2shm.obj.emp(par, U)

  par_sub <- l2shm.gd.proj(par, maxiter = 100, groups = 3)
  l2shm.obj.proj(par_sub, par)

  library(plotly)
  U <- rheat.sph.mix(10^3, par0)
  plot_ly() %>%
    add_trace(x = U[1, ], y = U[2, ], z = U[3, ], mode = "markers", marker=list(size=2)) %>%
    add_trace(x = 1.1*par$mu[1, ], y = 1.1*par$mu[2, ], z = 1.1*par$mu[3, ], mode = "markers", color=par$alpha)

  U_boot <- rheat.sph.mix(10^3, par)

  plot_ly() %>%
    add_trace(x = U[1, ], y = U[2, ], z = U[3, ], mode = "markers", marker=list(size=2)) %>%
    add_trace(x = U_boot[1, ], y = U_boot[2, ], z = U_boot[3, ], mode = "markers", marker=list(size=2)) %>%
    add_trace(x = 1.1*par$mu[1, ], y = 1.1*par$mu[2, ], z = 1.1*par$mu[3, ], mode = "markers", color=par$alpha)

  l2.heat.obj.emp(U, y)
}
