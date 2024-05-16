#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>


#include "l2shm_R.h"
#include "l2shm_random_R.h"

#include "l2shm.h"
#include "l2shm_random.h"


SEXP L2SHM_R(rand_heat_sph)(
    SEXP b_sxp,
    SEXP mu_sxp, SEXP T_sxp,
    SEXP res_sxp
)
{
    double * mu = REAL(mu_sxp);

    int n = Rf_xlength(mu_sxp);

    double T = Rf_asReal(T_sxp);

    int b = Rf_asInteger(b_sxp);
    int res = Rf_asInteger(res_sxp);

    SEXP out_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, b));
    double * out = REAL(out_sxp);

    for(size_t ind = 0; ind < b; ++ind)
        L2SHM(random_heat_sphere)(&out[ind*n], n, mu, T, res);

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(rand_heat_sph_mix)(
    SEXP b_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP res_sxp
)
{
    double * mu = REAL(mu_sxp);
    double * T = REAL(T_sxp);
    double * Alpha = REAL(Alpha_sxp);

    int n = Rf_nrows(mu_sxp);
    int k = Rf_ncols(mu_sxp);

    int b = Rf_asInteger(b_sxp);
    int res = Rf_asInteger(res_sxp);

    SEXP out_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, b));
    double * out = REAL(out_sxp);

    for(size_t ind = 0; ind < b; ++ind)
        L2SHM(random_heat_sphere_mixture)(&out[ind*n], n, k, mu, T, Alpha, res);

    UNPROTECT(1);
    return out_sxp;
}
