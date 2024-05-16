#ifndef L2SHM_GD_R_H_
#define L2SHM_GD_R_H_


#include "l2shm_R.h"


SEXP L2SHM_R(grad_desc_emp)
(
    SEXP U_sxp,
    SEXP t_min_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP max_iter_sxp, SEXP tol_sxp, SEXP terms_sxp
);

SEXP L2SHM_R(grad_desc_proj)
(
    SEXP mu0_sxp, SEXP T0_sxp, SEXP Alpha0_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP max_iter_sxp, SEXP tol_sxp, SEXP terms_sxp
);


#endif // L2SHM_GD_R_H_
