#ifndef L2SHM_MIXTURE_H_
#define L2SHM_MIXTURE_H_


#include "l2shm.h"


void L2SHM(copy_parameters)
(
    double *restrict mu, double *restrict T, double *restrict Alpha,
    size_t n, size_t k, size_t k0,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0
);

double L2SHM(objective_empirical)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double terms
);

double L2SHM(objective_projection)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
);

double L2SHM(norm2_squared)
(
    size_t n, size_t k,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double terms
);

double L2SHM(norm2_squared_difference)
(
    size_t n, size_t k1, size_t k2,
    double *restrict mu1, double *restrict T1, double *restrict Alpha1,
    double *restrict mu2, double *restrict T2, double *restrict Alpha2,
    double terms
);


#endif // L2SHM_MIXTURE_H_

