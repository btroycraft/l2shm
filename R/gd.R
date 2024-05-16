rand.par.list.emp <- function(U, groups, tmin)
{
  groups <- collapse.to.number(groups)
  tmin <- collapse.to.number(tmin)

  stopifnot(groups >= 1, tmin > 0)

    list(
        mu = U[, sample(ncol(U), groups)],
        t = rep(2*tmin, groups),
        alpha = rep(1/groups, groups)
    )
}

rand.par.list.proj <- function(par, groups)
{
  groups <- collapse.to.number(groups)

  stopifnot(groups >= 1)

  ind <- sample(ncol(par$mu), groups, prob = par$alpha)

    list(
        mu = par$mu[, ind],
        t = par$t[ind],
        alpha = rep(1/groups, groups)
    )
}

l2shm.gd.emp <- function(U, tmin, groups = 1L, par = rand.par.list.emp(U, groups, tmin), maxiter = 1L, tol = 10^-5, terms = 50L){
  invisible(.Call("_grad_desc_emp", PACKAGE="l2shm", U, tmin, par$mu, par$t, par$alpha, maxiter, tol, terms))
}

l2shm.gd.proj <- function(par0, groups = 1L, par = rand.par.list.proj(par0, groups), maxiter = 1L, tol = 10^-5, terms = 50L){
  invisible(.Call("_grad_desc_proj", PACKAGE="l2shm", par0$mu, par0$t, par0$alpha, par$mu, par$t, par$alpha, maxiter, tol, terms))
}
