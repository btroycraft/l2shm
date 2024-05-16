#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <math.h>


#include "l2shm_R.h"
#include "l2shm_gd_R.h"

#include "l2shm.h"
#include "l2shm_gd.h"


SEXP L2SHM_R(grad_desc_emp)
(
    SEXP U_sxp,
    SEXP T_min_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP max_iter_sxp, SEXP terms_sxp
)
{
    double * U = REAL(U_sxp);
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    int n = Rf_nrows(U_sxp);
    int p = Rf_ncols(U_sxp);
    int k = Rf_ncols(mu_sxp);

    double T_min = Rf_asReal(T_min_sxp);

    double terms = Rf_asReal(terms_sxp);
    int max_iter = Rf_asInteger(max_iter_sxp);

    double * alloc;
    {
        size_t length1 = 2*k + n;
        size_t length2 = 2*k + 2*k + k*(k+1)/2;
        alloc = (double *) R_alloc(length2 >= length1 ? length2 : length1, sizeof(double));
    }

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        T[ind] = exp(-T[ind]);

    L2SHM(gradient_descent_empirical)(n, k, p, mu, T, Alpha, U, exp(-T_min), alloc, terms, max_iter);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        T[ind] = -log(T[ind]);

    const char *out_names[] = {"mu", "t", "alpha", ""};

    SEXP out_sxp = PROTECT(Rf_mkNamed(VECSXP, out_names));

    SET_VECTOR_ELT(out_sxp, 0, mu_sxp);
    SET_VECTOR_ELT(out_sxp, 1, T_sxp);
    SET_VECTOR_ELT(out_sxp, 2, Alpha_sxp);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(grad_desc_proj)
(
    SEXP mu0_sxp, SEXP T0_sxp, SEXP Alpha0_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP max_iter_sxp, SEXP terms_sxp
)
{
    double * mu0 = REAL(mu0_sxp);
    double * T0 = REAL(T0_sxp);
    double * Alpha0 = REAL(Alpha0_sxp);
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    int n = Rf_nrows(mu0_sxp);
    int k0 = Rf_ncols(mu0_sxp);
    int k = Rf_ncols(mu_sxp);

    double terms = Rf_asReal(terms_sxp);
    int max_iter = Rf_asInteger(max_iter_sxp);

    double * alloc;
    {
        size_t length1 = 2*k + n;
        size_t length2 = 2*k + 2*k + k*(k+1)/2;
        alloc = (double *) R_alloc(length2 >= length1 ? length2 : length1, sizeof(double));
    }

    double * _T0 = (double *) R_alloc(k0, sizeof(double));

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k0; ++ind)
        _T0[ind] = exp(-T0[ind]);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        T[ind] = exp(-T[ind]);

    L2SHM(gradient_descent_projection)(n, k, k0, mu, T, Alpha, mu0, _T0, Alpha0, alloc, terms, max_iter);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        T[ind] = -log(T[ind]);

    const char *out_names[] = {"mu", "t", "alpha", ""};

    SEXP out_sxp = PROTECT(Rf_mkNamed(VECSXP, out_names));

    SET_VECTOR_ELT(out_sxp, 0, mu_sxp);
    SET_VECTOR_ELT(out_sxp, 1, T_sxp);
    SET_VECTOR_ELT(out_sxp, 2, Alpha_sxp);

    UNPROTECT(1);
    return out_sxp;
}
