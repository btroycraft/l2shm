#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <math.h>


#include "l2shm_R.h"
#include "l2shm_mixture_R.h"

#include "l2shm.h"
#include "l2shm_mixture.h"


SEXP L2SHM_R(copy_par)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP groups_sxp
)
{
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    int n = Rf_nrows(mu_sxp);
    int k = Rf_asInteger(groups_sxp);
    int k0 = Rf_ncols(mu_sxp);

    SEXP out_mu_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, k));
    SEXP out_T_sxp = PROTECT(Rf_allocVector(REALSXP, k));
    SEXP out_Alpha_sxp = PROTECT(Rf_allocVector(REALSXP, k));

    double * out_mu = REAL(out_mu_sxp);
    double * out_T = REAL(out_T_sxp);
    double * out_Alpha = REAL(out_Alpha_sxp);

    L2SHM(copy_parameters)(out_mu, out_T, out_Alpha, n, k, k0, mu, T, Alpha);

    const char *out_names[] = {"mu", "t", "alpha", ""};

    SEXP out_sxp = PROTECT(Rf_mkNamed(VECSXP, out_names));

    SET_VECTOR_ELT(out_sxp, 0, out_mu_sxp);
    SET_VECTOR_ELT(out_sxp, 1, out_T_sxp);
    SET_VECTOR_ELT(out_sxp, 2, out_Alpha_sxp);

    UNPROTECT(4);
    return out_sxp;
}

SEXP L2SHM_R(sort_par)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp
)
{
    const char *out_names[] = {"mu", "t", "alpha", ""};

    SEXP out_sxp = PROTECT(Rf_mkNamed(VECSXP, out_names));



    SET_VECTOR_ELT(out_sxp, 0, PROTECT(Rf_duplicate(mu_sxp)));
    SET_VECTOR_ELT(out_sxp, 1, PROTECT(Rf_duplicate(T_sxp)));
    SET_VECTOR_ELT(out_sxp, 2, PROTECT(Rf_duplicate(Alpha_sxp)));

    UNPROTECT(4);
    return out_sxp;
}

SEXP L2SHM_R(obj_emp)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP U_sxp,
    SEXP terms_sxp
)
{
    double * U = REAL(U_sxp);
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    size_t n = Rf_nrows(U_sxp);
    size_t p = Rf_ncols(U_sxp);
    size_t k = Rf_ncols(mu_sxp);

    double terms = Rf_asReal(terms_sxp);

    double * _T = (double *) R_alloc(k, sizeof(double));

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        _T[ind] = exp(-T[ind]);

    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));

    *REAL(out_sxp) = L2SHM(objective_empirical)(n, k, p, mu, _T, Alpha, U, terms);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(obj_proj)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP mu0_sxp, SEXP T0_sxp, SEXP Alpha0_sxp,
    SEXP terms_sxp
)
{
    double * mu0 = REAL(mu0_sxp);
    double * T0 = REAL(T0_sxp);
    double * Alpha0 = REAL(Alpha0_sxp);
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    size_t n = Rf_nrows(mu0_sxp);
    size_t k0 = Rf_ncols(mu0_sxp);
    size_t k = Rf_ncols(mu_sxp);

    double terms = Rf_asReal(terms_sxp);

    double * _T = (double *) R_alloc(k, sizeof(double));
    double * _T0 = (double *) R_alloc(k0, sizeof(double));

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        _T[ind] = exp(-T[ind]);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k0; ++ind)
        _T0[ind] = exp(-T0[ind]);

    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));

    *REAL(out_sxp) = L2SHM(objective_projection)(n, k, k0, mu, _T, Alpha, mu0, _T0, Alpha0, terms);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(nrm2_sq)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP terms_sxp
)
{
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    size_t n = Rf_nrows(mu_sxp);
    size_t k = Rf_ncols(mu_sxp);

    double terms = Rf_asReal(terms_sxp);

    double * _T = (double *) R_alloc(k, sizeof(double));

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k; ++ind)
        _T[ind] = exp(-T[ind]);

    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));

    *REAL(out_sxp) = L2SHM(norm2_squared)(n, k, mu, _T, Alpha, terms);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(nrm2_sq_diff)
(
    SEXP mu1_sxp, SEXP T1_sxp, SEXP Alpha1_sxp,
    SEXP mu2_sxp, SEXP T2_sxp, SEXP Alpha2_sxp,
    SEXP terms_sxp
)
{
    double * mu1 = REAL(mu1_sxp);
    double * T1 = REAL(T1_sxp);
    double * Alpha1 = REAL(Alpha1_sxp);
    double * mu2 = REAL(mu2_sxp);
    double * T2 = REAL(T2_sxp);
    double * Alpha2 = REAL(Alpha2_sxp);

    size_t n = Rf_nrows(mu1_sxp);
    size_t k1 = Rf_ncols(mu1_sxp);
    size_t k2 = Rf_ncols(mu2_sxp);

    double terms = Rf_asReal(terms_sxp);

    double * _T1 = (double *) R_alloc(k1, sizeof(double));
    double * _T2 = (double *) R_alloc(k2, sizeof(double));

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k1; ++ind)
        _T1[ind] = exp(-T1[ind]);

    #pragma GCC ivdep
    for(size_t ind = 0; ind < k2; ++ind)
        _T2[ind] = exp(-T2[ind]);

    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));

    *REAL(out_sxp) = L2SHM(norm2_squared_difference)(n, k1, k2, mu1, _T1, Alpha1, mu2, _T2, Alpha2, terms);

    UNPROTECT(1);
    return out_sxp;
}
