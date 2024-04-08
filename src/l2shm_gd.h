#ifndef L2SHM_GD_H_
#define L2SHM_GD_H_

#include "l2shm.h"

void L2SHM(grad_desc_emp)(size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double, double *restrict, double, size_t);
void L2SHM(grad_desc_proj)(size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double, double *restrict, double, size_t);

void L2SHM(stoch_grad_desc_emp)(size_t, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double, double *restrict, size_t *restrict, double, size_t);

#endif // L2SHM_GD_H_
