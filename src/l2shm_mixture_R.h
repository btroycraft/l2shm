#ifndef L2SHM_MIXTURE_R_H_
#define L2SHM_MIXTURE_R_H_


#include "l2shm_R.h"


SEXP L2SHM_R(copy_par)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP groups_sxp
);

SEXP L2SHM_R(obj_emp)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP U_sxp,
    SEXP terms_sxp
);

SEXP L2SHM_R(obj_proj)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP mu0_sxp, SEXP T0_sxp, SEXP Alpha0_sxp,
    SEXP terms_sxp
);

SEXP L2SHM_R(nrm2_sq)
(
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP terms_sxp
);

SEXP L2SHM_R(nrm2_sq_diff)
(
    SEXP mu1_sxp, SEXP T1_sxp, SEXP Alpha1_sxp,
    SEXP mu2_sxp, SEXP T2_sxp, SEXP Alpha2_sxp,
    SEXP terms_sxp
);


#endif // L2SHM_MIXTURE_R_H_

