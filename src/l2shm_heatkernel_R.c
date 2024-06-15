#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <math.h>

#include "l2shm_R.h"
#include "l2shm_heatkernel_R.h"

#include "l2shm.h"
#include "l2shm_heatkernel.h"


SEXP L2SHM_R(heat_kern)
(
    SEXP x_sxp, SEXP t_sxp,
    SEXP dim_sxp,
    SEXP terms_sxp
)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));

    *REAL(out_sxp) = L2SHM(heat_kernel)(Rf_asReal(x_sxp), exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(heat_kern_dx)
(
    SEXP x_sxp, SEXP t_sxp,
    SEXP dim_sxp,
    SEXP terms_sxp
)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));

    *REAL(out_sxp) = L2SHM(heat_kernel_dx)(Rf_asReal(x_sxp), exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(heat_kern_dt)
(
    SEXP x_sxp, SEXP t_sxp,
    SEXP dim_sxp,
    SEXP terms_sxp
)
{
    SEXP out_sxp = PROTECT(Rf_allocVector(REALSXP, 1));

    double T = exp(-Rf_asReal(t_sxp));

    *REAL(out_sxp) = -T * L2SHM(heat_kernel_dt)(Rf_asReal(x_sxp), T, Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(1);
    return out_sxp;
}

SEXP L2SHM_R(heat_kern_test)
(
    SEXP t_sxp, SEXP dim_sxp,
    SEXP terms_sxp
)
{
    SEXP out1_sxp = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP out2_sxp = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP out3_sxp = PROTECT(Rf_allocVector(REALSXP, 1));
    SEXP out4_sxp = PROTECT(Rf_allocVector(REALSXP, 2));

    SEXP out_sxp = PROTECT(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out_sxp, 0, out1_sxp);
    SET_VECTOR_ELT(out_sxp, 1, out2_sxp);
    SET_VECTOR_ELT(out_sxp, 2, out3_sxp);
    SET_VECTOR_ELT(out_sxp, 3, out4_sxp);

    double * out1 = REAL(out1_sxp);
    double * out2 = REAL(out2_sxp);
    double * out3 = REAL(out3_sxp);
    double * out4 = REAL(out4_sxp);

    out1[0] = L2SHM(heat_kernel)(1., exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));
    out1[1] = L2SHM(heat_kernel_dx)(1., exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));
    out1[2] = L2SHM(heat_kernel_dt)(1., exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    L2SHM(heat_kernel_combined_fdxdt)(&out2[0], &out2[1], &out2[2], 1., exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    out3[0] = L2SHM(heat_kernel_nox)(exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    L2SHM(heat_kernel_nox_combined_fdt)(&out4[0], &out4[1], exp(-Rf_asReal(t_sxp)), Rf_asInteger(dim_sxp), Rf_asReal(terms_sxp));

    UNPROTECT(5);
    return out_sxp;
}
