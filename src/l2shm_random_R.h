#ifndef L2SHM_RANDOM_R_H_
#define L2SHM_RANDOM_R_H_


#include "l2shm_R.h"


SEXP L2SHM_R(rand_heat_sph)(
    SEXP b_sxp,
    SEXP mu_sxp, SEXP T_sxp,
    SEXP res_sxp
);

SEXP L2SHM_R(rand_heat_sph_mix)(
    SEXP b_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP res_sxp
);


#endif // L2SHM_RANDOM_R_H_
