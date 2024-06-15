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
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double terms
);
static double L2SHM(group_objective)
(
    size_t index,
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double terms
);
static void L2SHM(Alpha_constants)
(
    double *restrict vec, double *restrict mat,
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T,
    double *restrict U,
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
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict t,
    double *restrict U,
    double *restrict mu_temp,
    double t_mult, double terms
);
static double L2SHM(update_Alpha)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict alpha,
    double *restrict U,
    double *restrict alpha_temp,
    double *restrict vec, double *restrict mat,
    double terms
);


void L2SHM(gradient_descent_empirical)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double t_mult,
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
        t[ind] = sqrt(t_mult/T[ind]-1.);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        alpha[ind] = sqrt(Alpha[ind]);

    double obj1;
    {
        for(size_t index = 0; index < k; ++index)
            L2SHM(update_group)(index, n, k, p, mu, T, Alpha, t, U, mu_temp, t_mult, terms);
        obj1 = L2SHM(update_Alpha)(n, k, p, mu, T, Alpha, alpha, U, alpha_temp, vec, mat, terms);
    }
    for(size_t iter = 1; iter < max_iter; ++iter){
        for(size_t index = 0; index < k; ++index)
            L2SHM(update_group)(index, n, k, p, mu, T, Alpha, t, U, mu_temp, t_mult, terms);
        double obj2 = L2SHM(update_Alpha)(n, k, p, mu, T, Alpha, alpha, U, alpha_temp, vec, mat, terms);
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
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double terms
)
{
    double *restrict mu_index = &mu[index*n];
    double T_index = T[index];

    double _obj, _T_grad;

// calculate term 1

    {
        double h, hx, ht;
        L2SHM(heat_kernel_combined_fdxdt)(&h, &hx, &ht, cblas_ddot(n, U, 1, mu_index, 1), T_index, n, terms);

        _obj = h;
        _T_grad = ht;
        #pragma GCC ivdep
        for(size_t ind = 0; ind < n; ++ind)
            mu_grad[ind] = hx*U[ind];
    }
    for(size_t ind = 1; ind < p; ++ind){
        double h, hx, ht;
        L2SHM(heat_kernel_combined_fdxdt)(&h, &hx, &ht, cblas_ddot(n, &U[ind*n], 1, mu_index, 1), T_index, n, terms);

        _obj += h;
        _T_grad += ht;
        cblas_daxpy(n, hx, &U[ind*n], 1, mu_grad, 1);
    }
    {
        double ip = 1./p;

        _obj *= ip;
        _T_grad *= ip;
        cblas_dscal(n, ip, mu_grad, 1);
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
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double terms
)
{
    double *restrict mu_index = &mu[index*n];
    double T_index = T[index];

// calculate term 1

    double obj = L2SHM(heat_kernel)(cblas_ddot(n, U, 1, mu_index, 1), T_index, n, terms);

    for(size_t ind = 1; ind < p; ++ind)
        obj += L2SHM(heat_kernel)(cblas_ddot(n, &U[ind*n], 1, mu_index, 1), T_index, n, terms);
    obj /= p;

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
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T,
    double *restrict U,
    double terms
)
{
    double ip = 1./p;
    for(size_t ind1 = 0; ind1 < k; ++ind1){
        double _T = T[ind1];
        double *restrict _mu = &mu[ind1*n];
        {
            double _obj = L2SHM(heat_kernel)(cblas_ddot(n, U, 1, _mu, 1), _T, n, terms);
            for(size_t ind2 = 1; ind2 < p; ++ind2)
                _obj += L2SHM(heat_kernel)(cblas_ddot(n, &U[ind2*n], 1, _mu, 1), _T, n, terms);
            vec[ind1] = _obj*ip;
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
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict t,
    double *restrict U,
    double *restrict mu_temp,
    double t_mult, double terms
)
{
    double *restrict mu_index = &mu[index*n];
    double t_index = t[index];

    double obj1, limit, t_temp;

    {
        double T_grad;
        L2SHM(group_combined_objective_gradient)(&obj1, mu_temp, &T_grad, index, n, k, p, mu, T, Alpha, U, terms);
        double t_grad = L2SHM(t_gradient_transform)(t_index, T_grad, t_mult);

        {
            double temp = cblas_dnrm2(n, mu_temp, 1);
            limit = .1*(temp*temp + t_grad*t_grad);
        }

        cblas_daxpy(n, 1., mu_index, 1, mu_temp, 1);
        cblas_dscal(n, 1./cblas_dnrm2(n, mu_temp, 1), mu_temp, 1);
        cblas_dswap(n, mu_index, 1, mu_temp, 1);

        t_temp = t_index;
        t_index = t_index + t_grad;
        T[index] = L2SHM(T_transform)(t_index, t_mult);
    }

    {
        double obj2 = L2SHM(group_objective)(index, n, k, p, mu, T, Alpha, U, terms);

        while(obj2 < obj1 + limit){

            limit *= .5;
            if(limit <= DBL_EPSILON){
                cblas_dcopy(n, mu_temp, 1, mu_index, 1);
                break;
            }

            cblas_daxpy(n, 1., mu_temp, 1, mu_index, 1);
            cblas_dscal(n, 1./cblas_dnrm2(n, mu_index, 1), mu_index, 1);

            t_index = .5*(t_temp + t_index);
            T[index] = L2SHM(T_transform)(t_index, t_mult);

            obj2 = L2SHM(group_objective)(index, n, k, p, mu, T, Alpha, U, terms);
        }
    }

    t[index] = t_index;

    return;
}

static double L2SHM(update_Alpha)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict alpha,
    double *restrict U,
    double *restrict alpha_temp,
    double *restrict vec, double *restrict mat,
    double terms
)
{
    L2SHM(Alpha_constants)(vec, mat, n, k, p, mu, T, U, terms);

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
