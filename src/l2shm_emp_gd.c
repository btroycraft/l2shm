#include <math.h>
#include <cblas.h>
#include <float.h>

#include <unistd.h>

#include <stdio.h>

#include "l2shm.h"
#include "l2shm_gd.h"
#include "l2shm_heatkernel.h"
#include "l2shm_trans.h"


static void L2SHM(group_combn)(double *, double *restrict, double *, size_t, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double);
static double L2SHM(group_obj)(size_t, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double);
static void L2SHM(Alpha_const)(double *restrict, double *restrict, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double);
static void L2SHM(Alpha_combn)(double *, double *restrict, size_t, double *restrict, double *restrict, double *restrict);
static double L2SHM(Alpha_obj)(size_t, double *restrict, double *restrict, double *restrict);
static double L2SHM(obj)(size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double);
static void L2SHM(update_group)(size_t, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double, double, double *);
static void L2SHM(update_Alpha)(size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double, double *);


void L2SHM(grad_desc_emp)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double t_min,
    double *restrict alloc,
    double terms, size_t max_iter
)
{
    double * t = alloc;
    double * alpha = t + k;

    double * mu_grad = alpha + k;
    double * mu_copy = mu_grad + n;

    double * alpha_grad = alpha + k;
    double * Alpha_grad = alpha_grad + k;
    double * alpha_copy = Alpha_grad + k;

    double * vec = alpha_copy + k;
    double * mat = vec + k;

    double t_mult = exp(-t_min);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        t[ind] = sqrt(exp(T[ind]-t_min)-1.);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind){
        T[ind] = L2SHM(T_trans)(t[ind], t_mult);
    }

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        alpha[ind] = sqrt(Alpha[ind]);

    double step = .5;
    for(size_t iter = 0; iter < max_iter; ++iter){
        printf("\nIteration %d: %e", (unsigned int) iter, L2SHM(obj)(n, k, p, mu, T, Alpha, U, terms));
        if(step <= DBL_EPSILON)
            step = .5;
        for(size_t index = 0; index < k; ++index)
            L2SHM(update_group)(index, n, k, p, mu, T, Alpha, t, U, mu_grad, mu_copy, t_mult, terms, &step);
        if(step <= DBL_EPSILON)
            step = .5;
        L2SHM(update_Alpha)(n, k, p, mu, T, Alpha, alpha, U, alpha_grad, Alpha_grad, alpha_copy, vec, mat, terms, &step);
    }

    for(size_t ind = 0; ind < k; ++ind)
        T[ind] = -log(T[ind]);

    return;
}

static void L2SHM(group_combn)
(
    double * obj, double *restrict mu_grad, double * T_grad,
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
        L2SHM(heatkernel_combn_dxdt)(&h, &hx, &ht, cblas_ddot(n, &U[0], 1, mu_index, 1), T_index, n, terms);

        _obj = h;
        _T_grad = ht;
        #pragma GCC ivdep
        for(size_t ind = 0; ind < n; ++ind)
            mu_grad[ind] = hx*U[ind];
    }
    for(size_t ind = 1; ind < p; ++ind){
        double h, hx, ht;
        L2SHM(heatkernel_combn_dxdt)(&h, &hx, &ht, cblas_ddot(n, &U[ind*n], 1, mu_index, 1), T_index, n, terms);

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
        double h, hx, ht;
        L2SHM(heatkernel_combn_dxdt)(&h, &hx, &ht, cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T_index*T[ind], n, terms);

        double _Alpha = Alpha[ind];
        _obj -= _Alpha*h;
        _T_grad -= _Alpha*T[ind]*ht;
        cblas_daxpy(n, -_Alpha*hx, &mu[ind*n], 1, mu_grad, 1);
    }
    for(size_t ind = index+1; ind < k; ++ind){
        double h, hx, ht;
        L2SHM(heatkernel_combn_dxdt)(&h, &hx, &ht, cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T_index*T[ind], n, terms);

        double _Alpha = Alpha[ind];
        _obj -= _Alpha*h;
        _T_grad -= _Alpha*T[ind]*ht;
        cblas_daxpy(n, -_Alpha*hx, &mu[ind*n], 1, mu_grad, 1);
    }

    // calculate term 3

    {
        double h, ht;
        L2SHM(heatkernel_combn_noxdt)(&h, &ht, T_index*T_index, n, terms);

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

static double L2SHM(group_obj)
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

    double obj = L2SHM(heatkernel)(cblas_ddot(n, &U[0], 1, mu_index, 1), T_index, n, terms);

    for(size_t ind = 1; ind < p; ++ind)
        obj += L2SHM(heatkernel)(cblas_ddot(n, &U[ind*n], 1, mu_index, 1), T_index, n, terms);
    obj /= p;

// calculate term 2

    for(size_t ind = 0; ind < index; ++ind)
        obj -= Alpha[ind]*L2SHM(heatkernel)(cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T_index*T[ind], n, terms);
    for(size_t ind = index+1; ind < k; ++ind)
        obj -= Alpha[ind]*L2SHM(heatkernel)(cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T_index*T[ind], n, terms);

    // calculate term 3

    obj -= .5*Alpha[index]*L2SHM(heatkernel_nox)(T_index*T_index, n, terms);

    return obj;
}

static void L2SHM(Alpha_const)
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
            double _obj = L2SHM(heatkernel)(cblas_ddot(n, &U[0], 1, _mu, 1), _T, n, terms);
            for(size_t ind2 = 1; ind2 < p; ++ind2)
                _obj += L2SHM(heatkernel)(cblas_ddot(n, &U[ind2*n], 1, _mu, 1), _T, n, terms);
            vec[ind1] = _obj*ip;
        }
        {
            double *restrict _mat = &mat[(ind1)*(ind1+1)/2];
            for(size_t ind2 = 0; ind2 < ind1; ++ind2)
                _mat[ind2] = L2SHM(heatkernel)(cblas_ddot(n, &mu[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
            _mat[ind1] = L2SHM(heatkernel_nox)(_T*_T, n, terms);
        }
    }

    return;
}

static void L2SHM(Alpha_combn)
(
    double * obj, double *restrict Alpha_grad,
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

static double L2SHM(Alpha_obj)
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

static double L2SHM(obj)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double terms
)
{
    double obj;
    {
        double *restrict _mu = &mu[0];
        double _T = T[0];
        double _obj = L2SHM(heatkernel)(cblas_ddot(n, &U[0], 1, _mu, 1), _T, n, terms);
        for(size_t ind = 1; ind < p; ++ind)
            _obj += L2SHM(heatkernel)(cblas_ddot(n, &U[ind*n], 1, _mu, 1), _T, n, terms);
        obj = Alpha[0]*_obj;
    }
    for(size_t ind1 = 1; ind1 < k; ++ind1){
        double *restrict _mu = &mu[ind1*n];
        double _T = T[ind1];
        double _obj = L2SHM(heatkernel)(cblas_ddot(n, &U[0], 1, _mu, 1), _T, n, terms);
        for(size_t ind2 = 1; ind2 < p; ++ind2)
            _obj += L2SHM(heatkernel)(cblas_ddot(n, &U[ind2*n], 1, _mu, 1), _T, n, terms);
        obj += Alpha[ind1]*_obj;
    }
    obj /= p;

    for(size_t ind1 = 0; ind1 < k; ++ind1){
        double *restrict _mu = &mu[ind1*n];
        double _Alpha = Alpha[ind1];
        double _T = T[ind1];
        double _obj = .5*_Alpha*L2SHM(heatkernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = ind1+1; ind2 < k; ++ind2)
            _obj += Alpha[ind2]*L2SHM(heatkernel)(cblas_ddot(n, &mu[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
        obj -= _Alpha*_obj;
    }

    return obj;
}

static void L2SHM(update_group)
(
    size_t index,
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict t,
    double *restrict U,
    double *restrict mu_grad, double *restrict mu_copy,
    double t_mult, double terms, double * step
)
{
    double t_index = t[index];
    double *restrict mu_index = &mu[index*n];

    double obj1, obj2, T_grad;

    L2SHM(group_combn)(&obj1, mu_grad, &T_grad, index, n, k, p, mu, T, Alpha, U, terms);

    double t_grad = L2SHM(t_grad_trans)(t_index, T_grad, t_mult);

    double mu_grad_norm = cblas_dnrm2(n, mu_grad, 1);
    double grad_norm2 = mu_grad_norm*mu_grad_norm + t_grad*t_grad;

    cblas_dcopy(n, mu_index, 1, mu_copy, 1);
    double t_copy = t_index;

    {
        double _step = *step;

        cblas_daxpy(n, _step, mu_grad, 1, mu_index, 1);
        cblas_dscal(n, 1./cblas_dnrm2(n, mu_index, 1), mu_index, 1);

        t_index = t_copy + _step*t_grad;
        T[index] = L2SHM(T_trans)(t_index, t_mult);

        obj2 = L2SHM(group_obj)(index, n, k, p, mu, T, Alpha, U, terms);

        if(obj2 >= obj1 + .5*_step*grad_norm2)
            _step *= 2.;
        else{
            do{
                if(_step <= DBL_EPSILON)
                    break;
                _step *= .5;

                cblas_dcopy(n, mu_copy, 1, mu_index, 1);

                cblas_daxpy(n, _step, mu_grad, 1, mu_index, 1);
                cblas_dscal(n, 1./cblas_dnrm2(n, mu_index, 1), mu_index, 1);

                t_index = t_copy + _step*t_grad;
                T[index] = L2SHM(T_trans)(t_index, t_mult);

                obj2 = L2SHM(group_obj)(index, n, k, p, mu, T, Alpha, U, terms);

            } while(obj2 < obj1 + .5*_step*grad_norm2);
        }

        *step = _step;
    }

    t[index] = t_index;

    return;
}

static void L2SHM(update_Alpha)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict alpha,
    double *restrict U,
    double *restrict alpha_grad, double *restrict Alpha_grad, double *restrict alpha_copy,
    double *restrict vec, double *restrict mat,
    double terms, double * step
)
{
    double obj1, obj2;

    L2SHM(Alpha_const)(vec, mat, n, k, p, mu, T, U, terms);

    L2SHM(Alpha_combn)(&obj1, Alpha_grad, k, Alpha, vec, mat);

    L2SHM(alpha_grad_trans)(alpha_grad, k, alpha, Alpha_grad);

    double grad_norm = cblas_dnrm2(k, alpha_grad, 1);
    double grad_norm2 = grad_norm*grad_norm;

    cblas_dcopy(k, alpha, 1, alpha_copy, 1);

    {
        double _step = *step;

        cblas_daxpy(k, _step, alpha_grad, 1, alpha, 1);
        cblas_dscal(k, 1./cblas_dnrm2(k, alpha, 1), alpha, 1);
        L2SHM(Alpha_trans)(Alpha, k, alpha);

        obj2 = L2SHM(Alpha_obj)(k, Alpha, vec, mat);

        if(obj2 >= obj1 + .5*_step*grad_norm2)
            _step *= 2.;
        else{
            do{
                if(_step <= DBL_EPSILON)
                    break;
                _step *= .5;

                cblas_dcopy(k, alpha_copy, 1, alpha, 1);

                cblas_daxpy(k, _step, alpha_grad, 1, alpha, 1);
                cblas_dscal(k, 1./cblas_dnrm2(k, alpha, 1), alpha, 1);
                L2SHM(Alpha_trans)(Alpha, k, alpha);

                obj2 = L2SHM(Alpha_obj)(k, Alpha, vec, mat);

            } while(obj2 < obj1 + .5*_step*grad_norm2);
        }

        *step = _step;
    }

    return;
}
