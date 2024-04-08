rand.par.list <- function(groups, n, tmin)
{
  groups <- collapse.to.number(groups)
  n <- collapse.to.number(n)
  tmin <- collapse.to.number(tmin)

  stopifnot(groups >= 1, n >= 1, tmin >= 0)

    list(
        mu = replicate(groups, {
              ._ <- rnorm(n)
              ._ / sqrt(sum(._^2))
            }),
        t = rep(2*tmin, groups),
        alpha = rep(1/groups, groups)
    )
}

heatkern.proj.gd <- function(par, tmin, groups, start = rand.par.list(groups, nrow(U), tmin), maxiter = 100L, terms = 50L){

  .Call("_heatkern_proj_gd", PACKAGE="l2shm", par, tmin, start, maxiter, terms)
}

heatkern.emp.gd <- function(U, tmin, groups = 1L, start = rand.par.list(groups, nrow(U), tmin), maxiter = 100L, terms = 50L){
  .Call("_heatkern_emp_gd", PACKAGE="l2shm", U, tmin, start, maxiter, terms)
}

heatkern.proj.sgd <- function(par, tmin, groups, batch, start = rand.par.list(groups, nrow(U), tmin), maxiter = 100L, terms = 50L){
  .Call("_heatkern_proj_sgd", PACKAGE="l2shm", par, tmin, batch, start, maxiter, terms)
}

heatkern.emp.sgd <- function(U, tmin, groups, batch, start = rand.par.list(groups, nrow(U), tmin), maxiter = 100L, terms = 50L){
  .Call("_heatkern_emp_sgd", PACKAGE="l2shm", U, tmin, batch, start, maxiter, terms)
}

heatkern <- function(x, t, dim, terms = 100L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heatkern", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heatkern.test <- function(x, t, dim, terms = 100L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heatkern_test", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heatkern.dx <- function(x, t, dim, terms = 100L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heatkern_dx", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heatkern.dx.test <- function(x, t, dim, terms = 100L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heatkern_dx_test", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heatkern.dt <- function(x, t, dim, terms = 100L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heatkern_dt", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heatkern.dt.test <- function(x, t, dim, terms = 100L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heatkern_dt_test", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heatkern.combn <- function(x, t, dim, terms = 100L){
    ._ <- sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heatkern_combn", PACKAGE="l2shm", .x, .t, dim, terms)
        }, c("f"=0, "dx"=0, "dt"=0))
    }, simplify = "array")
    list(f = ._["f", , ], dx = ._["dx", , ], dt = ._["dt", , ])
}
