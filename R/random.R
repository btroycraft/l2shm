rheat.sph <- function(n, par, res=30L){
    .Call("_rand_heat_sph", PACKAGE="l2shm", n, par$mu, par$t, res)
}

rheat.sph.mix <- function(n, par, res=30L){
    .Call("_rand_heat_sph_mix", PACKAGE="l2shm", n, par$mu, par$t, par$alpha, res)
}
