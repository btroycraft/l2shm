#ifndef L2SHM_RANDOM_H_
#define L2SHM_RANDOM_H_


#include "l2shm.h"


void L2SHM(random_batch)
(
    double *restrict U_batch,
    size_t n, size_t p, size_t b,
    double *restrict U,
    size_t *restrict indices
);

void L2SHM(random_heat_sphere)
(
    double *restrict out,
    size_t n,
    double *restrict mu, double T,
    size_t res
);

void L2SHM(random_heat_sphere_mixture)
(
    double *restrict out,
    size_t n, size_t k,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    size_t res
);


#endif // L2SHM_RANDOM_H_
