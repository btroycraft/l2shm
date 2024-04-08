#include <math.h>
#include <cblas.h>


#include "l2shm.h"
#include "l2shm_gd.h"
#include "l2shm_heatkernel.h"
#include "l2shm_trans.h"


static void L2SHM(group_combn)(double *, double *restrict, double *, size_t, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double);
static double L2SHM(group_obj)(size_t, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double);
static void L2SHM(alpha_const)(double *restrict, double *restrict, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double);
static void L2SHM(alpha_combn)(double *, double *restrict, size_t, double *restrict, double *restrict, double *restrict);
static double L2SHM(alpha_obj)(size_t, double *restrict, double *restrict, double *restrict);
static void L2SHM(update_group)(size_t, size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double, double, double *, size_t *, size_t);
static void L2SHM(update_alpha)(size_t, size_t, size_t, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double *restrict, double, double *, size_t *, size_t);


void L2SHM(grad_desc_proj)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double t_min,
    double *restrict alloc,
    double terms, size_t max_iter
)
{
    double * t = &alloc[0];
    double * alpha = &alloc[k];
    double * mu_grad = &alloc[2*k];
    double * mu_copy = &alloc[2*k + n];
    double * alpha_grad = &alloc[2*k];
    double * Alpha_grad = &alloc[3*k];
    double * alpha_copy = &alloc[4*k];
    double * vec = &alloc[5*k];
    double * mat = &alloc[6*k];

    double t_mult = exp(-t_min);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        t[ind] = sqrt(exp(T[ind]-t_min)-1.);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        T[ind] = L2SHM(T_trans)(t[ind], t_mult);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        alpha[ind] = sqrt(Alpha[ind]);

    double step = 1.;
    for(size_t iter = 0; iter < max_iter; ++iter){
        for(size_t index = 0; index < k; ++index)
            L2SHM(update_group)(index, n, k, k0, mu, T, Alpha, t, mu0, T0, Alpha0, mu_grad, mu_copy, t_mult, terms, &step, &iter, max_iter);
        L2SHM(update_alpha)(n, k, k0, mu, T, Alpha, alpha, mu0, T0, Alpha0, Alpha_grad, alpha_grad, alpha_copy, vec, mat, terms, &step, &iter, max_iter);
    }

    for(size_t ind = 0; ind < k; ++ind)
        T[ind] = -log(T[ind]);

    return;
}


static void L2SHM(group_combn)
(
    double * obj, double *restrict mu_grad, double * T_grad,
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
        double _Alpha0 = Alpha0[0];
        double h, hx, ht;
        L2SHM(heatkernel_combn_dxdt)(&h, &hx, &ht, cblas_ddot(n, &mu0[0], 1, mu_index, 1), T_index*T0[0], n, terms);

        _obj = _Alpha0*h;
        _T_grad = _Alpha0*ht;
        {
            double temp = _Alpha0*hx;
            #pragma GCC ivdep
            for(size_t ind = 0; ind < k; ++ind)
                mu_grad[ind] = temp*mu0[ind];
        }
    }
    for(size_t ind = 1; ind < k0; ++ind){

        double _Alpha0 = Alpha0[ind];
        double h, hx, ht;
        L2SHM(heatkernel_combn_dxdt)(&h, &hx, &ht, cblas_ddot(n, &mu0[ind*n], 1, mu_index, 1), T_index*T0[ind], n, terms);

        _obj += _Alpha0*h;
        _T_grad += _Alpha0*ht;
        cblas_daxpy(n, _Alpha0*hx, &mu0[ind*n], 1, mu_grad, 1);
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
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
)
{
    double *restrict mu_index = &mu[index*n];
    double T_index = T[index];

// calculate term 1

    double obj = Alpha0[0]*L2SHM(heatkernel)(cblas_ddot(n, &mu0[0], 1, mu_index, 1), T_index*T0[0], n, terms);

    for(size_t ind = 1; ind < k0; ++ind)
        obj += Alpha0[ind]*L2SHM(heatkernel)(cblas_ddot(n, &mu0[ind*n], 1, mu_index, 1), T_index*T0[ind], n, terms);

// calculate term 2

    for(size_t ind = 0; ind < index; ++ind)
        obj -= Alpha[ind]*L2SHM(heatkernel)(cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T_index*T[ind], n, terms);
    for(size_t ind = index+1; ind < k; ++ind)
        obj -= Alpha[ind]*L2SHM(heatkernel)(cblas_ddot(n, &mu[ind*n], 1, mu_index, 1), T_index*T[ind], n, terms);

    // calculate term 3

    obj -= .5*Alpha[index]*L2SHM(heatkernel_nox)(T_index*T_index, n, terms);

    return obj;
}

static void L2SHM(alpha_const)
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
            double _obj = Alpha0[0]*L2SHM(heatkernel)(cblas_ddot(n, &mu0[0], 1, _mu, 1), _T*T[0], n, terms);
            for(size_t ind2 = 1; ind2 < k0; ++ind2)
                _obj += Alpha0[ind2]*L2SHM(heatkernel)(cblas_ddot(n, &mu0[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
            vec[ind1] = _obj;
        }
        {
            double *restrict _mat = &mat[(ind1)*(ind1+1)/2];
            _mat[ind1] = L2SHM(heatkernel_nox)(_T*_T, n, terms);
            for(size_t ind2 = 1; ind2 <= ind1; ++ind2)
                _mat[ind2] = L2SHM(heatkernel)(cblas_ddot(n, &mu[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
        }
    }

    return;
}

static void L2SHM(alpha_combn)
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

static double L2SHM(alpha_obj)
(
    size_t k,
    double *restrict Alpha,
    double *restrict vec, double *restrict mat
)
{
    double obj = cblas_ddot(k, Alpha, 1, vec, 1);
    obj -= .5*Alpha[0]*mat[0];
    for(size_t ind = 1; ind < k; ++ind){
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
    double *restrict mu_grad, double *restrict mu_copy,
    double t_mult, double terms, double * step, size_t * iter, size_t max_iter
)
{
    double t_index = t[index];
    double *restrict mu_index = &mu[index*n];

    double obj1=0., obj2=0., T_grad=0.;

    L2SHM(group_combn)(&obj1, mu_grad, &T_grad, index, n, k, k0, mu, T, Alpha, mu0, T0, Alpha0, terms);

    double t_grad = L2SHM(t_grad_trans)(t_index, T_grad, t_mult);

    double mu_grad_norm = cblas_dnrm2(n, mu_grad, 1);
    double grad_norm2 = mu_grad_norm*mu_grad_norm + t_grad*t_grad;

    cblas_dcopy(n, mu_index, 1, mu_copy, 1);
    double t_copy = t_index;

    size_t _iter = *iter;
    double _step = *step;
    for(; _iter < max_iter; ++_iter){

        cblas_daxpy(n, _step, mu_grad, 1, mu_index, 1);
        cblas_dscal(n, 1./cblas_dnrm2(n, mu_index, 1), mu_index, 1);

        t_index = t_copy + _step*t_grad;
        T[index] = L2SHM(T_trans)(t_index, t_mult);

        obj2 = L2SHM(group_obj)(index, n, k, k0, mu, T, Alpha, mu0, T0, Alpha0, terms);

        if(obj2 >= obj1)
            break;

        _step *= .9;

        cblas_dcopy(n, mu_copy, 1, mu_index, 1);
    }
    *iter = _iter;
    *step = _step;

    double diff = obj2-obj1;

    if(diff < _step*grad_norm2){
        double mult = _step*_step*grad_norm2/(2.*(_step*grad_norm2-diff));

        t_index = t_copy + mult*t_grad;
        T[index] = L2SHM(T_trans)(t_index, t_mult);

        cblas_dcopy(n, mu_copy, 1, mu_index, 1);
        cblas_daxpy(n, mult, mu_grad, 1, mu_index, 1);
        cblas_dscal(n, 1./cblas_dnrm2(n, mu_index, 1), mu_index, 1);
    }

    return;
}

static void L2SHM(update_alpha)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha, double *restrict alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double *restrict Alpha_grad, double *restrict alpha_grad, double *restrict alpha_copy,
    double *restrict vec, double *restrict mat,
    double terms, double * step, size_t * iter, size_t max_iter
)
{
    double obj1 = 0., obj2 = 0. ;

    L2SHM(alpha_const)(vec, mat, n, k, k0, mu, T, mu0, T0, Alpha0, terms);

    L2SHM(alpha_combn)(&obj1, Alpha_grad, k, Alpha, vec, mat);
    L2SHM(alpha_grad_trans)(alpha_grad, k, alpha, Alpha_grad);

    double grad_norm = cblas_dnrm2(k, alpha_grad, 1);
    double grad_norm2 = grad_norm*grad_norm;

    cblas_dcopy(k, alpha, 1, alpha_copy, 1);

    size_t _iter = *iter;
    double _step = *step;
    for(; _iter < max_iter; ++_iter){

        cblas_daxpy(k, _step, alpha_grad, 1, alpha, 1);
        cblas_dscal(k, 1./cblas_dnrm2(k, alpha, 1), alpha, 1);
        L2SHM(Alpha_trans)(Alpha, k, alpha);

        obj2 = L2SHM(alpha_obj)(k, Alpha, vec, mat);

        if(obj2 >= obj1)
            break;

        _step *= .9;

        cblas_dcopy(k, alpha_copy, 1, alpha, 1);
    }
    *iter = _iter;
    *step = _step;

    double diff = obj2-obj1;

    if(diff < _step*grad_norm2){
        double mult = _step*_step*grad_norm2/(2.*(_step*grad_norm2-diff));

        cblas_dcopy(k, alpha_copy, 1, alpha, 1);
        cblas_daxpy(k, mult, alpha_grad, 1, alpha, 1);
        cblas_dscal(k, 1./cblas_dnrm2(n, alpha, 1), alpha, 1);
    }

    return;
}
