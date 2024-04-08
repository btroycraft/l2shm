u.scores <- function(X){
  X <- as.matrix(X)
  n <- nrow(X)
  Y <- sweep(X[-n, ], 2, (colSums(X) + sqrt(n)*X[n, ])/(n+sqrt(n)))
  sweep(Y, 2, sqrt(colSums(Y^2)), '/')
}
