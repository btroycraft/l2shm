#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>


#include "l2shm_R.h"
#include "l2shm_random_R.h"

#include "l2shm.h"
#include "l2shm_random.h"


SEXP L2SHM_R(rand_unif_sph)
(
    SEXP n_sxp,
    SEXP dim_sxp
)
{
    int n = Rf_asInteger(n_sxp);
    int dim = Rf_asInteger(dim_sxp);

    SEXP out_sxp = PROTECT(Rf_allocMatrix(REALSXP, dim, n));
    double * out = REAL(out_sxp);

    for(size_t ind = 0; ind < n; ++ind)
        L2SHM(random_uniform_sphere)(&out[ind*dim], dim);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(rand_heat_sph)
(
    SEXP n_sxp,
    SEXP mu_sxp, SEXP T_sxp,
    SEXP res_sxp
)
{
    double * mu = REAL(mu_sxp);

    int dim = Rf_xlength(mu_sxp);

    double T = Rf_asReal(T_sxp);

    int n = Rf_asInteger(n_sxp);
    int res = Rf_asInteger(res_sxp);

    SEXP out_sxp = PROTECT(Rf_allocMatrix(REALSXP, dim, n));
    double * out = REAL(out_sxp);

    for(size_t ind = 0; ind < n; ++ind)
        L2SHM(random_heat_sphere)(&out[ind*dim], dim, mu, T, res);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(rand_heat_sph_mix)
(
    SEXP n_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP res_sxp
)
{
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    int dim = Rf_nrows(mu_sxp);
    int k = Rf_ncols(mu_sxp);

    int n = Rf_asInteger(n_sxp);
    int res = Rf_asInteger(res_sxp);

    SEXP out_sxp = PROTECT(Rf_allocMatrix(REALSXP, dim, n));
    double * out = REAL(out_sxp);

    for(size_t ind = 0; ind < n; ++ind)
        L2SHM(random_heat_sphere_mixture)(&out[ind*dim], dim, k, mu, T, Alpha, res);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(rand_combn)
(
    SEXP n_sxp,
    SEXP k_sxp,
    SEXP max_sxp
)
{
    int n = Rf_asInteger(n_sxp);
    int k = Rf_asInteger(k_sxp);
    int max = Rf_asInteger(max_sxp);

    size_t * indices = (size_t *) R_alloc(k, sizeof(size_t));

    SEXP out_sxp = PROTECT(Rf_allocMatrix(INTSXP, k, n));
    int * out = INTEGER(out_sxp);

    for(size_t ind1 = 0; ind1 < n; ++ind1){
        L2SHM(random_combination(indices, k, max-1));
        int *restrict _out = &out[ind1*k];
        #pragma GCC ivdep
        for(size_t ind2 = 0; ind2 < k; ++ind2)
            _out[ind2] = (int) indices[ind2] + 1;
    }

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(rand_int)
(
    SEXP n_sxp,
    SEXP max_sxp
)
{
    int n = Rf_asInteger(n_sxp);
    int max = Rf_asInteger(max_sxp);

    SEXP out_sxp = PROTECT(Rf_allocVector(INTSXP, n));
    int * out = INTEGER(out_sxp);

    for(size_t ind = 0; ind < n; ++ind)
        out[ind] = (int) L2SHM(random_integer(max-1) + 1);

    UNPROTECT(1);
    return out_sxp;
}
