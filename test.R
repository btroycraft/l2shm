install.packages("./", repos=NULL, type="source")
library(l2shm)

U <- {
  data <- read.csv("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv", header=TRUE)
  data <- subset(data, location == 'Westchester Lagoon', c(-location, -sample))
  u.scores(data)
}

start <- rand.par.list(10, nrow(U), .02)

y <- heatkern.emp.gd(U, tmin = .02, groups = 10, maxiter = 20)
y <- heatkern.emp.gd(U, .02, 10, maxiter = 1000, start=y); y

y_max <- heatkern.emp.sgd(U, .02, 4, 10, maxiter = 10)

for(._ in 1:20){
  y <- heatkern.emp.sgd(U, .02, 4, 10, maxiter = 20)
  y <- heatkern.emp.sgd(U, .02, maxiter = 200, start=y, batch = 10)
  y <- heatkern.emp.gd(U, .02, maxiter = 100, start=y)
  y <- heatkern.emp.sgd(U, .02, maxiter = 200, start=y, batch = 10)
}


plot(seq(-1, 1, .01), heatkern(seq(-1, 1, .01), exp(-.02), 23, 27))
plot(seq(-1, 1, .01), heatkern.dx(seq(-1, 1, .01), .95, 30, 50))

t <- head(seq(0, 1, .01), -1)


{._ <- runif(4); ._/sum(._)}

{
  alpha <- local((._ <- runif(4))/sum(._))
  mu <- replicate(length(alpha), (._ <- rnorm(3))/sqrt(sum(._^2)))

  X <- matrix(rnorm(nrow(mu)*100, sd=.2), nrow(mu), 100)
  U <- local({
    ._ <- X+mu[, sample(ncol(mu), ncol(X), prob = alpha, replace = TRUE)]
    sweep(._, 2, sqrt(colSums(._^2)), '/')
  })

  y <- heatkern.emp.sgd(U, .0001, 100, 20, maxiter = 10)
  y <- heatkern.emp.sgd(U, .0001, 100, 20, maxiter = , start = y)
  y <- heatkern.emp.gd(U, .0001, 100, maxiter = 10^3, start = y)

  library(plotly)
  plot_ly() %>%
    add_trace(x = U[1, ], y = U[2, ], z = U[3, ], mode = "markers") %>%
    add_trace(x = mu[1, ], y = mu[2, ], z = mu[3, ], mode = "markers") %>%
    add_trace(x = y$mu[1, ], y = y$mu[2, ], z = y$mu[3, ], mode = "markers", color=y$alpha)

}
