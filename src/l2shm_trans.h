#ifndef L2SHM_TRANS_H_
#define L2SHM_TRANS_H_


#include "l2shm.h"

double L2SHM(T_trans)(double, double);
double L2SHM(t_grad_trans)(double, double, double);
void L2SHM(Alpha_trans)(double *restrict, size_t, double *restrict);
void L2SHM(alpha_grad_trans)(double *restrict, size_t, double *restrict, double *restrict);


#endif // L2SHM_TRANS_H_
