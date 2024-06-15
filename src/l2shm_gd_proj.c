#include <math.h>
#include <cblas.h>
#include <float.h>
#include <unistd.h>


#include "l2shm.h"
#include "l2shm_gd.h"
#include "l2shm_heatkernel.h"
#include "l2shm_mixture.h"
#include "l2shm_trans.h"


static void L2SHM(group_combined_objective_gradient)
(
    double *restrict obj, double *restrict mu_grad, double *restrict T_grad,
    size_t index,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
);
static double L2SHM(group_objective)
(
    size_t index,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
);
static void L2SHM(Alpha_constants)
(
    double *restrict vec, double *restrict mat,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
);
static void L2SHM(Alpha_combined_objective_gradient)
(
    double *restrict obj, double *restrict Alpha_grad,
    size_t k,
    double *restrict Alpha,
    double *restrict vec, double *restrict mat
);
static double L2SHM(Alpha_objective)
(
    size_t k,
    double *restrict Alpha,
    double *restrict vec, double *restrict mat
);
static void L2SHM(update_group)
(
    size_t index,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict t,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double *restrict mu_temp,
    double terms
);
static double L2SHM(update_Alpha)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double *restrict alpha_temp,
    double *restrict vec, double *restrict mat,
    double terms
);


void L2SHM(gradient_descent_projection)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double * alloc,
    double terms, size_t max_iter, double tol
)
{
    double * t = alloc;
    double * alpha = t + k;

    double * mu_temp = alpha + k;

    double * alpha_temp = alpha + k;

    double * vec = alpha_temp + k;
    double * mat = vec + k;

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        t[ind] = sqrt(1./T[ind]-1.);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        alpha[ind] = sqrt(Alpha[ind]);

    double obj1;
    {
        for(size_t index = 0; index < k; ++index)
            L2SHM(update_group)(index, n, k, k0, mu, T, Alpha, t, mu0, T0, Alpha0, mu_temp, terms);
        obj1 = L2SHM(update_Alpha)(n, k, k0, mu, T, Alpha, alpha, mu0, T0, Alpha0, alpha_temp, vec, mat, terms);
    }
    for(size_t iter = 1; iter < max_iter; ++iter){
        for(size_t index = 0; index < k; ++index)
            L2SHM(update_group)(index, n, k, k0, mu, T, Alpha, t, mu0, T0, Alpha0, mu_temp, terms);
        double obj2 = L2SHM(update_Alpha)(n, k, k0, mu, T, Alpha, alpha, mu0, T0, Alpha0, alpha_temp, vec, mat, terms);
        if(obj2-obj1 <= tol)
            break;
        obj1 = obj2;
    }

    return;
}

static void L2SHM(group_combined_objective_gradient)
(
    double *restrict obj, double *restrict mu_grad, double *restrict T_grad,
    size_t index,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
)
{
    double *restrict mu_index = &mu[index*n];
    double T_index = T[index];

    double _obj, _T_grad;

// calculate term 1

    {
        double _T = T0[0];
        double h, hx, ht;
        L2SHM(heat_kernel_combined_fdxdt)(&h, &hx, &ht, cblas_ddot(n, mu0, 1, mu_index, 1), _T*T_index, n, terms);

        double _Alpha = Alpha0[0];
        _obj = _Alpha*h;
        _T_grad = _Alpha*_T*ht;
        {
            double temp = _Alpha*hx;
            #pragma GCC ivdep
            for(size_t ind = 0; ind < n; ++ind)
                mu_grad[ind] = temp*mu0[ind];
        }
    }
    for(size_t ind = 1; ind < k0; ++ind){
        double _T = T0[ind];
        double h, hx, ht;
        L2SHM(heat_kernel_combined_fdxdt)(&h, &hx, &ht, cblas_ddot(n, &mu0[ind*n], 1, mu_index, 1), _T*T_index, n, terms);

        double _Alpha = Alpha0[ind];
        _obj += _Alpha*h;
        _T_grad += _Alpha*_T*ht;
        cblas_daxpy(n, _Alpha*hx, &mu0[ind*n], 1, mu_grad, 1);
    }

// calculate term 2


    for(size_t ind = 0; ind < index; ++ind){
        double _T = T[ind];
        double h, hx, ht;
        L2SHM(heat_kernel_combined_fdxdt)(&h, &hx, &ht, cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), _T*T_index, n, terms);

        double _Alpha = Alpha[ind];
        _obj -= _Alpha*h;
        _T_grad -= _Alpha*_T*ht;
        cblas_daxpy(n, -_Alpha*hx, &mu[ind*n], 1, mu_grad, 1);
    }
    for(size_t ind = index+1; ind < k; ++ind){
        double _T = T[ind];
        double h, hx, ht;
        L2SHM(heat_kernel_combined_fdxdt)(&h, &hx, &ht, cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), _T*T_index, n, terms);

        double _Alpha = Alpha[ind];
        _obj -= _Alpha*h;
        _T_grad -= _Alpha*_T*ht;
        cblas_daxpy(n, -_Alpha*hx, &mu[ind*n], 1, mu_grad, 1);
    }

    // calculate term 3

    {
        double h, ht;
        L2SHM(heat_kernel_nox_combined_fdt)(&h, &ht, T_index*T_index, n, terms);

        double _Alpha = Alpha[index];
        _obj -= .5*_Alpha*h;
        _T_grad -= _Alpha*T_index*ht;
    }

// Orthogonalize mu gradient

    cblas_daxpy(n, -cblas_ddot(n, mu_grad, 1, mu_index, 1), mu_index, 1, mu_grad, 1);

// Output

    *obj = _obj;
    *T_grad = _T_grad;

    return;
}

static double L2SHM(group_objective)
(
    size_t index,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
)
{
    double *restrict mu_index = &mu[index*n];
    double T_index = T[index];

// calculate term 1

    double obj = Alpha0[0]*L2SHM(heat_kernel)(cblas_ddot(n, mu0, 1, mu_index, 1), T0[0]*T_index, n, terms);

    for(size_t ind = 1; ind < k0; ++ind)
        obj += Alpha0[ind]*L2SHM(heat_kernel)(cblas_ddot(n, &mu0[ind*n], 1, mu_index, 1), T0[ind]*T_index, n, terms);

// calculate term 2

    for(size_t ind = 0; ind < index; ++ind)
        obj -= Alpha[ind]*L2SHM(heat_kernel)(cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T[ind]*T_index, n, terms);
    for(size_t ind = index+1; ind < k; ++ind)
        obj -= Alpha[ind]*L2SHM(heat_kernel)(cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T[ind]*T_index, n, terms);

    // calculate term 3

    obj -= .5*Alpha[index]*L2SHM(heat_kernel_nox)(T_index*T_index, n, terms);

    return obj;
}

static void L2SHM(Alpha_constants)
(
    double *restrict vec, double *restrict mat,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
)
{
    for(size_t ind1 = 0; ind1 < k; ++ind1){
        double _T = T[ind1];
        double *restrict _mu = &mu[ind1*n];
        {
            double _obj = Alpha0[0]*L2SHM(heat_kernel)(cblas_ddot(n, mu0, 1, _mu, 1), _T*T0[0], n, terms);
            for(size_t ind2 = 1; ind2 < k0; ++ind2){
                _obj += Alpha0[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu0[ind2*n], 1, _mu, 1), _T*T0[ind2], n, terms);
            }

            vec[ind1] = _obj;
        }
        {
            double *restrict _mat = &mat[ind1*(ind1+1)/2];
            for(size_t ind2 = 0; ind2 < ind1; ++ind2)
                _mat[ind2] = L2SHM(heat_kernel)(cblas_ddot(n, &mu[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
            _mat[ind1] = L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        }
    }

    return;
}

static void L2SHM(Alpha_combined_objective_gradient)
(
    double *restrict obj, double *restrict Alpha_grad,
    size_t k,
    double *restrict Alpha,
    double *restrict vec, double *restrict mat
)
{
    cblas_dcopy(k, vec, 1, Alpha_grad, 1);
    cblas_dspmv(CblasColMajor, CblasUpper, k, -1., mat, Alpha, 1, 1., Alpha_grad, 1);

    *obj = .5*(cblas_ddot(k, Alpha, 1, Alpha_grad, 1) + cblas_ddot(k, Alpha, 1, vec, 1));

    return;
}

static double L2SHM(Alpha_objective)
(
    size_t k,
    double *restrict Alpha,
    double *restrict vec, double *restrict mat
)
{
    double obj = cblas_ddot(k, Alpha, 1, vec, 1);
    for(size_t ind = 0; ind < k; ++ind){
        double *restrict _mat = &mat[ind*(ind+1)/2];
        double _Alpha = Alpha[ind];
        obj -= _Alpha*(cblas_ddot(ind, Alpha, 1, _mat, 1)+.5*_Alpha*_mat[ind]);
    }

    return obj;
}

static void L2SHM(update_group)
(
    size_t index,
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict t,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double *restrict mu_temp,
    double terms
)
{
    double *restrict mu_index = &mu[index*n];
    double t_index = t[index];

    double obj1, limit, t_temp;

    {
        double T_grad;
        L2SHM(group_combined_objective_gradient)(&obj1, mu_temp, &T_grad, index, n, k, k0, mu, T, Alpha, mu0, T0, Alpha0, terms);
        double t_grad = L2SHM(t_gradient_transform_nomult)(t_index, T_grad);

        {
            double temp = cblas_dnrm2(n, mu_temp, 1);
            limit = .1*(temp*temp + t_grad*t_grad);
        }

        cblas_daxpy(n, 1., mu_index, 1, mu_temp, 1);
        cblas_dscal(n, 1./cblas_dnrm2(n, mu_temp, 1), mu_temp, 1);
        cblas_dswap(n, mu_index, 1, mu_temp, 1);

        t_temp = t_index;
        t_index = t_index + t_grad;
        T[index] = L2SHM(T_transform_nomult)(t_index);
    }

    {
        double obj2 = L2SHM(group_objective)(index, n, k, k0, mu, T, Alpha, mu0, T0, Alpha0, terms);

        while(obj2 < obj1 + limit){

            limit *= .5;
            if(limit <= DBL_EPSILON){
                cblas_dcopy(n, mu_temp, 1, mu_index, 1);
                break;
            }

            cblas_daxpy(n, 1., mu_temp, 1, mu_index, 1);
            cblas_dscal(n, 1./cblas_dnrm2(n, mu_index, 1), mu_index, 1);

            t_index = .5*(t_temp + t_index);
            T[index] = L2SHM(T_transform_nomult)(t_index);

            obj2 = L2SHM(group_objective)(index, n, k, k0, mu, T, Alpha, mu0, T0, Alpha0, terms);
        }
    }

    t[index] = t_index;

    return;
}

static double L2SHM(update_Alpha)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double *restrict alpha_temp,
    double *restrict vec, double *restrict mat,
    double terms
)
{
    L2SHM(Alpha_constants)(vec, mat, n, k, k0, mu, T, mu0, T0, Alpha0, terms);

    double obj1, limit;

    {
        L2SHM(Alpha_combined_objective_gradient)(&obj1, alpha_temp, k, Alpha, vec, mat);
        L2SHM(alpha_gradient_transform)(alpha_temp, k, alpha);

        {
            double temp = cblas_dnrm2(k, alpha_temp, 1);
            limit = .1*temp*temp;
        }

        cblas_daxpy(k, 1., alpha, 1, alpha_temp, 1);
        cblas_dscal(k, 1./cblas_dnrm2(k, alpha_temp, 1), alpha_temp, 1);
        cblas_dswap(k, alpha, 1, alpha_temp, 1);

        L2SHM(Alpha_transform)(Alpha, k, alpha);
    }

    double obj2;
    {
        obj2 = L2SHM(Alpha_objective)(k, Alpha, vec, mat);

        while(obj2 < obj1 + limit){

            limit *= .5;
            if(limit <= DBL_EPSILON){
                cblas_dcopy(k, alpha_temp, 1, alpha, 1);
                L2SHM(Alpha_transform)(Alpha, k, alpha);
                obj2 = obj1;
                break;
            }

            cblas_daxpy(k, 1., alpha_temp, 1, alpha, 1);
            cblas_dscal(k, 1./cblas_dnrm2(k, alpha, 1), alpha, 1);

            L2SHM(Alpha_transform)(Alpha, k, alpha);

            obj2 = L2SHM(Alpha_objective)(k, Alpha, vec, mat);
        }
    }

    return obj2;
}
