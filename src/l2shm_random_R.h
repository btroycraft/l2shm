#ifndef L2SHM_RANDOM_R_H_
#define L2SHM_RANDOM_R_H_

#include <R.h>
#include <Rinternals.h>

#include "l2shm_R.h"


SEXP L2SHM_R(rand_unif_sph)
(
    SEXP n_sxp,
    SEXP dim_sxp
);

SEXP L2SHM_R(rand_heat_sph)
(
    SEXP n_sxp,
    SEXP mu_sxp, SEXP T_sxp,
    SEXP res_sxp
);

SEXP L2SHM_R(rand_heat_sph_mix)
(
    SEXP n_sxp,
    SEXP mu_sxp, SEXP T_sxp, SEXP Alpha_sxp,
    SEXP res_sxp
);

SEXP L2SHM_R(rand_combn)
(
    SEXP n_sxp,
    SEXP k_sxp,
    SEXP max_sxp
);

SEXP L2SHM_R(rand_int)
(
    SEXP n_sxp,
    SEXP max_sxp
);


#endif // L2SHM_RANDOM_R_H_
