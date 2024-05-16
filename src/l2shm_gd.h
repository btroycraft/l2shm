#ifndef L2SHM_GD_H_
#define L2SHM_GD_H_


#include "l2shm.h"


void L2SHM(gradient_descent_empirical)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double t_mult,
    double *restrict alloc,
    double terms, size_t max_iter, double tol
);

void L2SHM(gradient_descent_projection)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double *restrict alloc,
    double terms, size_t max_iter, double tol
);

#endif // L2SHM_GD_H_
