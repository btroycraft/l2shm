collapse.to.element <- function(x){
  while(is.list(x)){
    stopifnot(length(x) == 1)
    x <- x[[1]]
  }
  return(unname(x))
}

collapse.to.number <- function(x){
  x <- collapse.to.element(x)
  stopifnot(is.atomic(x) && is.numeric(x))
  stopifnot(length(x) == 1)
  return(x[1])
}
