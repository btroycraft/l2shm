l2shm.obj.emp <- function(par, U, terms = 30L){
    .Call("_obj_emp", PACKAGE="l2shm", par$mu, par$t, par$alpha, U, terms)
}

l2shm.obj.proj <- function(par, par0, terms = 30L){
    .Call("_obj_proj", PACKAGE="l2shm", par$mu, par$t, par$alpha, par0$mu, par0$t, par0$alpha, terms)
}

l2shm.nrm2.sq <- function(par, par0, terms = 30L){
    .Call("_nrm2_sq", PACKAGE="l2shm", par$mu, par$t, par$alpha, terms)
}

l2shm.nrm2.sq.diff <- function(par1, par2, terms = 30L){
    .Call("_nrm2_sq_diff", PACKAGE="l2shm", par1$mu, par1$t, par1$alpha, par2$mu, par2$t, par2$alpha, terms)
}
