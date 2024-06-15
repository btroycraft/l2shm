rint <- function(n, max){
    .Call("_rand_int", PACKAGE="l2shm", n, max)
}

rcombn <- function(n, k, max){
    .Call("_rand_combn", PACKAGE="l2shm", n, k, max)
}

runif.sph <- function(n, dim){
    .Call("_rand_unif_sph", PACKAGE="l2shm", n, dim)
}

rheat.sph <- function(n, par, res=30L){
    .Call("_rand_heat_sph", PACKAGE="l2shm", n, par$mu, par$t, res)
}

rheat.sph.mix <- function(n, par, res=30L){
    .Call("_rand_heat_sph_mix", PACKAGE="l2shm", n, par$mu, par$t, par$alpha, res)
}
