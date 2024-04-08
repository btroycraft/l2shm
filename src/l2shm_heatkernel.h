#ifndef L2SHM_L2SHMEL_H_
#define L2SHM_L2SHMEL_H_

#include "l2shm.h"


double L2SHM(heatkernel)(double, double, size_t, double);
double L2SHM(heatkernel_nox)(double, size_t, double);
double L2SHM(heatkernel_dx)(double, double, size_t, double);
double L2SHM(heatkernel_dxx)(double, double, size_t, double);
double L2SHM(heatkernel_dt)(double, double, size_t, double);
void L2SHM(heatkernel_combn_dxdt)(double *, double *, double *, double, double, size_t, double);
void L2SHM(heatkernel_combn_noxdt)(double *, double *, double, size_t, double);
double L2SHM(heatkernel_test)(double, double, size_t, double);
double L2SHM(heatkernel_dx_test)(double, double, size_t, double);
double L2SHM(heatkernel_dxx_test)(double, double, size_t, double);
double L2SHM(heatkernel_dt_test)(double, double, size_t, double);


#endif // L2SHM_L2SHMEL_H_
