heat.kern <- function(x, t, dim, terms = 30L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heat_kern", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heat.kern.dx <- function(x, t, dim, terms = 30L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heat_kern_dx", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}

heat.kern.dt <- function(x, t, dim, terms = 30L){
    sapply(t, function(.t){
        vapply(x, function(.x){
            .Call("_heat_kern_dt", PACKAGE="l2shm", .x, .t, dim, terms)
        }, numeric(1))
    })
}
