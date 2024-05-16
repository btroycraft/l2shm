#ifndef L2SHM_HEATKERNEL_R_H_
#define L2SHM_HEATKERNEL_R_H_


#include "l2shm_R.h"


SEXP L2SHM_R(heat_kern)
(
    SEXP x_sxp, SEXP t_sxp,
    SEXP dim_sxp,
    SEXP terms_sxp
);

SEXP L2SHM_R(heat_kern_dx)
(
    SEXP x_sxp, SEXP t_sxp,
    SEXP dim_sxp,
    SEXP terms_sxp
);

SEXP L2SHM_R(heat_kern_dt)
(
    SEXP x_sxp, SEXP t_sxp,
    SEXP dim_sxp,
    SEXP terms_sxp
);

SEXP L2SHM_R(heat_kern_test)
(
    SEXP t_sxp, SEXP dim_sxp,
    SEXP terms_sxp
);


#endif // L2SHM_HEATKERNEL_R_H_
