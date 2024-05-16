#ifndef L2SHM_TRANS_H_
#define L2SHM_TRANS_H_


#include "l2shm.h"


double L2SHM(T_transform)
(
    double t,
    double t_mult
);

double L2SHM(T_transform_nomult)
(
    double t
);

double L2SHM(t_gradient_transform)
(
    double t, double T_grad,
    double t_mult
);

double L2SHM(t_gradient_transform_nomult)
(
    double t, double T_grad
);

void L2SHM(Alpha_transform)
(
    double *restrict Alpha,
    size_t k,
    double *restrict alpha
);

void L2SHM(alpha_gradient_transform)
(
    double *restrict alpha_grad,
    size_t k,
    double *restrict alpha
);


#endif // L2SHM_TRANS_H_
