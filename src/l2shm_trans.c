#include <cblas.h>


#include "l2shm.h"
#include "l2shm_trans.h"


double L2SHM(T_transform)
(
    double t,
    double t_mult
)
{
    return t_mult / (1.+t*t);
}

double L2SHM(T_transform_nomult)
(
    double t
)
{
    return 1. / (1.+t*t);
}

double L2SHM(t_gradient_transform)
(
    double t, double T_grad,
    double t_mult
)
{
    double temp = 1.+t*t;
    return -2.*t_mult*t*T_grad / (temp*temp);
}

double L2SHM(t_gradient_transform_nomult)
(
    double t, double T_grad
)
{
    double temp = 1.+t*t;
    return -2.*t*T_grad / (temp*temp);
}

void L2SHM(Alpha_transform)
(
    double *restrict Alpha,
    size_t k,
    double *restrict alpha
)
{
    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind){
        double _alpha = alpha[ind];
        Alpha[ind] = _alpha*_alpha;
    }

    return;
}

void L2SHM(alpha_gradient_transform)
(
    double *restrict alpha_grad,
    size_t k,
    double *restrict alpha
)
{
    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        alpha_grad[ind] *= 2.*alpha[ind];

    cblas_daxpy(k, -cblas_ddot(k, alpha_grad, 1, alpha, 1), alpha, 1, alpha_grad, 1);

    return;
}
