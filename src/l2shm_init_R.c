#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


#include "l2shm_init_R.h"


static const R_CallMethodDef callMethods[] = {
    {"_grad_desc_emp", (DL_FUNC) &L2SHM_R(grad_desc_emp), 8},
    {"_grad_desc_proj", (DL_FUNC) &L2SHM_R(grad_desc_proj), 9},
    {"_obj_emp", (DL_FUNC) &L2SHM_R(obj_emp), 5},
    {"_obj_proj", (DL_FUNC) &L2SHM_R(obj_proj), 7},
    {"_nrm2_sq", (DL_FUNC) &L2SHM_R(nrm2_sq), 4},
    {"_nrm2_sq_diff", (DL_FUNC) &L2SHM_R(nrm2_sq_diff), 7},
    {"_rand_heat_sph", (DL_FUNC) &L2SHM_R(rand_heat_sph), 4},
    {"_rand_heat_sph_mix", (DL_FUNC) &L2SHM_R(rand_heat_sph_mix), 5},
    {"_heat_kern", (DL_FUNC) &L2SHM_R(heat_kern), 4},
    {"_heat_kern_dx", (DL_FUNC) &L2SHM_R(heat_kern_dx), 4},
    {"_heat_kern_dt", (DL_FUNC) &L2SHM_R(heat_kern_dt), 4},
    {NULL, NULL, 0}
};

void R_init_l2shm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
