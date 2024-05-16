#ifndef L2SHM_HEATKERNEL_H_
#define L2SHM_HEATKERNEL_H_


#include "l2shm.h"


double L2SHM(heat_kernel)
(
    double x, double t,
    size_t dim,
    double terms
);
double L2SHM(heat_kernel_nox)
(
    double t,
    size_t dim,
    double terms
);
double L2SHM(heat_kernel_dx)
(
    double x, double t,
    size_t dim,
    double terms
);
double L2SHM(heat_kernel_dxx)
(
    double x, double t,
    size_t dim,
    double terms
);
double L2SHM(heat_kernel_dt)
(
    double x, double t,
    size_t dim,
    double terms
);
void L2SHM(heat_kernel_combined_dxdt)
(
    double *restrict out, double *restrict out_dx, double *restrict out_dt,
    double x, double t,
    size_t dim,
    double terms
);
void L2SHM(heat_kernel_nox_combined_dt)
(
    double *restrict out, double *restrict out_dt,
    double t,
    size_t dim,
    double terms
);
double L2SHM(heat_kernel_test)
(
    double x, double t,
    size_t dim,
    double terms
);


#endif // L2SHM_HEATKERNEL_H_
