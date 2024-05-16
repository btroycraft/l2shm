#include <math.h>
#include <cblas.h>

#include "l2shm.h"
#include "l2shm_mixture.h"
#include "l2shm_heatkernel.h"


double L2SHM(objective_empirical)
(
    size_t n, size_t k, size_t p,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict U,
    double terms
)
{
    double obj;
    {
        double _T = T[0];
        double _obj = L2SHM(heat_kernel)(cblas_ddot(n, U, 1, mu, 1), _T, n, terms);
        for(size_t ind = 1; ind < p; ++ind)
            _obj += L2SHM(heat_kernel)(cblas_ddot(n, &U[ind*n], 1, mu, 1), _T, n, terms);
        obj = Alpha[0]*_obj;
    }
    for(size_t ind1 = 1; ind1 < k; ++ind1){
        double *restrict _mu = &mu[ind1*n];
        double _T = T[ind1];
        double _obj = L2SHM(heat_kernel)(cblas_ddot(n, U, 1, _mu, 1), _T, n, terms);
        for(size_t ind2 = 1; ind2 < p; ++ind2)
            _obj += L2SHM(heat_kernel)(cblas_ddot(n, &U[ind2*n], 1, _mu, 1), _T, n, terms);
        obj += Alpha[ind1]*_obj;
    }
    obj /= p;

    for(size_t ind1 = 0; ind1 < k; ++ind1){
        double *restrict _mu = &mu[ind1*n];
        double _Alpha = Alpha[ind1];
        double _T = T[ind1];
        double _obj = .5*_Alpha*L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = ind1+1; ind2 < k; ++ind2)
            _obj += Alpha[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
        obj -= _Alpha*_obj;
    }

    return obj;
}

double L2SHM(objective_projection)
(
    size_t n, size_t k, size_t k0,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double *restrict mu0, double *restrict T0, double *restrict Alpha0,
    double terms
)
{
    double obj;
    {
        double _T = T[0];
        double _obj = Alpha0[0]*L2SHM(heat_kernel)(cblas_ddot(n, mu0, 1, mu, 1), _T*T0[0], n, terms);
        for(size_t ind = 1; ind < k0; ++ind)
            _obj += Alpha0[ind]*L2SHM(heat_kernel)(cblas_ddot(n, &mu0[ind*n], 1, mu, 1), _T*T0[ind], n, terms);
        obj = Alpha[0]*_obj;
    }
    for(size_t ind1 = 1; ind1 < k; ++ind1){
        double *restrict _mu = &mu[ind1*n];
        double _T = T[ind1];
        double _obj = Alpha0[0]*L2SHM(heat_kernel)(cblas_ddot(n, mu0, 1, _mu, 1), _T*T0[0], n, terms);
        for(size_t ind2 = 1; ind2 < k0; ++ind2)
            _obj += Alpha0[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu0[ind2*n], 1, _mu, 1), _T*T0[ind2], n, terms);
        obj += Alpha[ind1]*_obj;
    }

    for(size_t ind1 = 0; ind1 < k; ++ind1){
        double *restrict _mu = &mu[ind1*n];
        double _Alpha = Alpha[ind1];
        double _T = T[ind1];
        double _obj = .5*_Alpha*L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = ind1+1; ind2 < k; ++ind2)
            _obj += Alpha[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
        obj -= _Alpha*_obj;
    }

    return obj;
}

double L2SHM(norm2_squared)
(
    size_t n, size_t k,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    double terms
)
{
    double obj;

    {
        double _Alpha = Alpha[0];
        double _T = T[0];
        double _obj = .5*_Alpha*L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = 1; ind2 < k; ++ind2)
            _obj += Alpha[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu[ind2*n], 1, mu, 1), _T*T[ind2], n, terms);
        obj = _Alpha*_obj;
    }
    for(size_t ind1 = 1; ind1 < k; ++ind1){
        double *restrict _mu = &mu[ind1*n];
        double _Alpha = Alpha[ind1];
        double _T = T[ind1];
        double _obj = .5*_Alpha*L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = ind1+1; ind2 < k; ++ind2)
            _obj += Alpha[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu[ind2*n], 1, _mu, 1), _T*T[ind2], n, terms);
        obj += _Alpha*_obj;
    }

    return obj;
}

double L2SHM(norm2_squared_difference)
(
    size_t n, size_t k1, size_t k2,
    double *restrict mu1, double *restrict T1, double *restrict Alpha1,
    double *restrict mu2, double *restrict T2, double *restrict Alpha2,
    double terms
)
{
    double obj;

    {
        double _Alpha = Alpha1[0];
        double _T = T1[0];
        double _obj = .5*_Alpha*L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = 1; ind2 < k1; ++ind2)
            _obj += Alpha1[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu1[ind2*n], 1, mu1, 1), _T*T1[ind2], n, terms);
        obj = _Alpha*_obj;
    }
    for(size_t ind1 = 1; ind1 < k1; ++ind1){
        double *restrict _mu = &mu1[ind1*n];
        double _Alpha = Alpha1[ind1];
        double _T = T1[ind1];
        double _obj = .5*_Alpha*L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = ind1+1; ind2 < k1; ++ind2)
            _obj += Alpha1[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu1[ind2*n], 1, _mu, 1), _T*T1[ind2], n, terms);
        obj += _Alpha*_obj;
    }

    for(size_t ind1 = 0; ind1 < k2; ++ind1){
        double *restrict _mu = &mu2[ind1*n];
        double _Alpha = Alpha2[ind1];
        double _T = T2[ind1];
        double _obj = .5*_Alpha*L2SHM(heat_kernel_nox)(_T*_T, n, terms);
        for(size_t ind2 = ind1+1; ind2 < k2; ++ind2)
            _obj += Alpha2[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu2[ind2*n], 1, _mu, 1), _T*T2[ind2], n, terms);
        obj += _Alpha*_obj;
    }

    for(size_t ind1 = 0; ind1 < k1; ++ind1){
        double *restrict _mu = &mu1[ind1*n];
        double _T = T1[ind1];
        double _obj = Alpha2[0]*L2SHM(heat_kernel)(cblas_ddot(n, mu2, 1, _mu, 1), _T*T2[0], n, terms);
        for(size_t ind2 = 1; ind2 < k2; ++ind2)
            _obj += Alpha2[ind2]*L2SHM(heat_kernel)(cblas_ddot(n, &mu2[ind2*n], 1, _mu, 1), _T*T2[ind2], n, terms);
        obj -= Alpha1[ind1]*_obj;
    }

    return obj;
}
