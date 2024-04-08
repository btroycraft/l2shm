#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <stdlib.h>

#include "l2shm_R.h"
#include "l2shm_random_R.h"

#include "l2shm.h"
#include "l2shm_heatkernel.h"
#include "l2shm_gd.h"

SEXP heatkern(SEXP x_sxp, SEXP t_sxp, SEXP dim_sxp, SEXP terms_sxp)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));
    *REAL(out_sxp) = L2SHM(heatkernel)(Rf_asReal(x_sxp), Rf_asReal(t_sxp), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}


SEXP heatkern_test(SEXP x_sxp, SEXP t_sxp, SEXP dim_sxp, SEXP terms_sxp)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));
    *REAL(out_sxp) = L2SHM(heatkernel_test)(Rf_asReal(x_sxp), Rf_asReal(t_sxp), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP heatkern_dx(SEXP x_sxp, SEXP t_sxp, SEXP dim_sxp, SEXP terms_sxp)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));
    *REAL(out_sxp) = L2SHM(heatkernel_dx)(Rf_asReal(x_sxp), Rf_asReal(t_sxp), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP heatkern_dx_test(SEXP x_sxp, SEXP t_sxp, SEXP dim_sxp, SEXP terms_sxp)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));
    *REAL(out_sxp) = L2SHM(heatkernel_dx_test)(Rf_asReal(x_sxp), Rf_asReal(t_sxp), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP heatkern_dt(SEXP x_sxp, SEXP t_sxp, SEXP dim_sxp, SEXP terms_sxp)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));
    *REAL(out_sxp) = L2SHM(heatkernel_dt)(Rf_asReal(x_sxp), Rf_asReal(t_sxp), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP heatkern_dt_test(SEXP x_sxp, SEXP t_sxp, SEXP dim_sxp, SEXP terms_sxp)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));
    *REAL(out_sxp) = L2SHM(heatkernel_dt_test)(Rf_asReal(x_sxp), Rf_asReal(t_sxp), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP heatkern_combn(SEXP x_sxp, SEXP t_sxp, SEXP dim_sxp, SEXP terms_sxp)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 3));
    double *restrict out = REAL(out_sxp);

    L2SHM(heatkernel_combn_dxdt)(&out[0], &out[1], &out[2], Rf_asReal(x_sxp), Rf_asReal(t_sxp), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP heatkern_emp_gd(SEXP U_sxp, SEXP t_min_sxp, SEXP start_sxp, SEXP max_iter_sxp, SEXP terms_sxp)
{
    double *restrict U = REAL(U_sxp);
    size_t n = Rf_nrows(U_sxp);
    size_t p = Rf_ncols(U_sxp);

    double t_min = Rf_asReal(t_min_sxp);

    double terms = Rf_asReal(terms_sxp);
    size_t max_iter = Rf_asInteger(max_iter_sxp);

    SEXP mu_sxp = VECTOR_ELT(start_sxp, 0);
    SEXP T_sxp = VECTOR_ELT(start_sxp, 1);
    SEXP Alpha_sxp = VECTOR_ELT(start_sxp, 2);

    double *restrict mu = REAL(mu_sxp);
    double *restrict T = REAL(T_sxp);
    double *restrict Alpha = REAL(Alpha_sxp);
    size_t k = Rf_ncols(mu_sxp);

    double * alloc;
    {
        size_t length1 = 2*k + 2*n;
        size_t length2 = 6*k + k*(k+1)/2;
        alloc = (double *) R_alloc(length2 >= length1 ? length2 : length1, sizeof(double));
    }

    L2SHM(grad_desc_emp)(n, k, p, mu, T, Alpha, U, t_min, alloc, terms, max_iter);

    return start_sxp;
}

SEXP heatkern_emp_sgd(SEXP U_sxp, SEXP t_min_sxp, SEXP batch_sxp, SEXP start_sxp, SEXP max_iter_sxp, SEXP terms_sxp)
{
    double *restrict U = REAL(U_sxp);
    size_t n = Rf_nrows(U_sxp);
    size_t p = Rf_ncols(U_sxp);

    double t_min = Rf_asReal(t_min_sxp);

    size_t b = Rf_asInteger(batch_sxp);
    double terms = Rf_asReal(terms_sxp);
    size_t max_iter = Rf_asInteger(max_iter_sxp);

    SEXP mu_sxp = VECTOR_ELT(start_sxp, 0);
    SEXP T_sxp = VECTOR_ELT(start_sxp, 1);
    SEXP Alpha_sxp = VECTOR_ELT(start_sxp, 2);

    double *restrict mu = REAL(mu_sxp);
    double *restrict T = REAL(T_sxp);
    double *restrict Alpha = REAL(Alpha_sxp);
    size_t k = Rf_ncols(mu_sxp);

    double * alloc1;
    {
        size_t length1 = 2*k + n + b*n;
        size_t length2 = 5*k + k*(k+1)/2 + b*n;
        alloc1 = (double *) R_alloc(length2 >= length1 ? length2 : length1, sizeof(double));
    }

    size_t * alloc2 = (size_t *) R_alloc(b, sizeof(size_t));

    L2SHM(stoch_grad_desc_emp)(n, k, p, b, mu, T, Alpha, U, t_min, alloc1, alloc2, terms, max_iter);

    return start_sxp;
}

SEXP heatkern_proj_gd(SEXP par_sxp, SEXP t_min_sxp, SEXP start_sxp, SEXP max_iter_sxp, SEXP terms_sxp)
{
    SEXP mu0_sxp = VECTOR_ELT(par_sxp, 0);
    SEXP T0_sxp = VECTOR_ELT(par_sxp, 1);
    SEXP Alpha0_sxp = VECTOR_ELT(par_sxp, 2);

    double *restrict mu0 = REAL(mu0_sxp);
    double *restrict T0 = REAL(T0_sxp);
    double *restrict Alpha0 = REAL(Alpha0_sxp);
    size_t n = Rf_nrows(mu0_sxp);
    size_t k0 = Rf_ncols(mu0_sxp);

    double t_min = Rf_asReal(t_min_sxp);

    double terms = Rf_asReal(terms_sxp);
    size_t max_iter = Rf_asInteger(max_iter_sxp);

    SEXP mu_sxp = VECTOR_ELT(start_sxp, 0);
    SEXP T_sxp = VECTOR_ELT(start_sxp, 1);
    SEXP Alpha_sxp = VECTOR_ELT(start_sxp, 2);

    double *restrict mu = REAL(mu_sxp);
    double *restrict T = REAL(T_sxp);
    double *restrict Alpha = REAL(Alpha_sxp);
    size_t k = Rf_ncols(mu_sxp);

    double * alloc;
    {
        size_t length1 = 2*k+n;
        size_t length2 = k*(k+13)/2;
        alloc = (double *) R_alloc(length2 >= length1 ? length2 : length1, sizeof(double));
    }

    L2SHM(grad_desc_proj)(n, k, k0, mu, T, Alpha, mu0, T0, Alpha0, t_min, alloc, terms, max_iter);

    return start_sxp;
}

#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"_heatkern", (DL_FUNC) &heatkern, 4},
    {"_heatkern_test", (DL_FUNC) &heatkern_test, 4},
    {"_heatkern_dx", (DL_FUNC) &heatkern_dx, 4},
    {"_heatkern_dx_test", (DL_FUNC) &heatkern_dx_test, 4},
    {"_heatkern_dt", (DL_FUNC) &heatkern_dt, 4},
    {"_heatkern_dt_test", (DL_FUNC) &heatkern_dt_test, 4},
    {"_heatkern_combn", (DL_FUNC) &heatkern_combn, 4},
    {"_heatkern_emp_gd", (DL_FUNC) &heatkern_emp_gd, 5},
    {"_heatkern_emp_sgd", (DL_FUNC) &heatkern_emp_sgd, 6},
    {"_heatkern_proj_gd", (DL_FUNC) &heatkern_proj_gd, 5},
    {NULL, NULL, 0}
};

void R_init_l2shm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
