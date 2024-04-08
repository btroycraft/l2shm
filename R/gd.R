reduce.to.element <- function(x){
  while(is.list(x)){
    stopifnot(length(x) == 1)
    x <- x[[1]]
  }
  return(x)
}

reduce.to.number <- function(x){
  while(is.list(x)){
    stopifnot(length(x) == 1)
    x <- x[[1]]
  }
  stopifnot(is.atomic(x) && is.numeric(x))
  stopifnot(length(x) == 1)
  return(x[1])
}

rand.par.list <- function(groups, n, tmin)
{
    {
        stopifnot(is.vector(groups), is.vector(n), is.vector(tmin))

        if(is.list(groups)) groups <- unlist(groups)
        if(is.list(n)) n <- unlist(n)
        if(is.list(tmin)) tmin <- unlist(tmin)

        stopifnot(
            is.numeric(groups) && is.atomic(groups),
            is.numeric(n) && is.atomic(n),
            is.numeric(t) && is.atomic(t)
        )
        stopifnot(length(groups) == 1, length(n) == 1, length(tmin) == 1)

        groups <- groups[1]
        n <- n[1]
        tmin <- tmin[1]

        stopifnot(groups >= 1, n >= 1, tmin >= 0)
    }

    list(
        mu = {
            ._ <- replicate()
            ._ <- matrix(rnorm(groups*n), n, groups)
            sweep(._, 2, sqrt(colSums(._^2)), '/')
        },
        t = rep(2*tmin, groups),
        alpha = rep(1/groups, groups)
    )
}

heatkern.proj.gd <- function(par, tmin, groups, start = rand.par.list(groups, nrow(U), tmin), maxiter = 100L, terms = 50L){

  .Call("_heatkern_proj_gd", PACKAGE="l2shm", par, tmin, start, maxiter, terms)
}

heatkern.emp.gd <- function(U, tmin, groups = 1L, start = rand.par.list(groups, nrow(U), tmin), maxiter = 100L, terms = 50L){
    {
        stopifnot(is.vector(tmin), is.vector(groups), is.vector(maxiter), is.vector(terms))

        if(is.list(tmin)) tmin <- unlist(tmin)
        if(is.list(groups)) groups <- unlist(groups)
        if(is.list(maxiter)) maxiter <- unlist(maxiter)
        if(is.list(terms)) terms <- unlist(terms)

        stopifnot(
            is.numeric(tmin) && is.atomic(tmin),
            is.numeric(groups) && is.atomic(groups),
            is.numeric(maxiter) && is.atomic(maxiter),
            is.numeric(terms) && is.atomic(terms)
        )
        stopifnot(length(tmin) == 1, length(groups) == 1, length(maxiter) == 1, length(terms) == 1)

        tmin <- tmin[1]
        groups <- groups[1]
        maxiter <- maxiter[1]
        terms <- terms[1]

        stopifnot(tmin >= 0, groups >= 1, maxiter >= 1, n >= 1)
    }

    {
        stopifnot(is.vector(U))

        if(is.list(U)) U <- unlist(U)

        stopifnot(is.numeric(U) && is.atomic(U))

        if(!is.matrix(U)) U <- as.matrix(U)

        {
            ._ <- apply(unname(U), 2, function(._) sum(._^2))

        }
        stopifnot()
    }


    stopifnot(
        is.matrix(U) && is.numeric(U),
        is.numeric(U)

    )



    stopifnot(
        is.numeric(maxiter) && length(maxiter) == 1


    )
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
