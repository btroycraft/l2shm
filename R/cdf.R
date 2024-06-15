cdf.lin <- function(x){
    x_sort <- sort(x)
    function(y){
        .Call("_cdf_lin", PACKAGE = "l2shm", sort(y), x_sort)[rank(y, ties.method = "first")]
    }
}

icdf <- function(x){
    x_sort <- sort(x)
    function(y){
        out <- sort(y)
        .Call("_cdf_inv_lin", PACKAGE = "l2shm", out, x_sort)

    }
}
