#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <stdint.h>


#include "l2shm_R.h"
#include "l2shm_cdf_R.h"

#include "l2shm.h"
#include "l2shm_cdf.h"


SEXP L2SHM_R(cdf_lin)
(
    SEXP y_sxp,
    SEXP x_sort_sxp
)
{
    double * y = REAL(y_sxp);
    double * x_sort = REAL(x_sort_sxp);

    size_t length_y = Rf_xlength(y_sxp);
    size_t length_x_sort = Rf_xlength(x_sort_sxp);

    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, length_y));
    double * out = REAL(out_sxp);

    L2SHM(cdf_linear)(out, length_y, length_x_sort, y, x_sort);

    UNPROTECT(1);
    return out_sxp;
}
