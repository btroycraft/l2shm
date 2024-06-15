#ifndef L2SHM_CDF_H_
#define L2SHM_CDF_H_


#include "l2shm.h"


void L2SHM(cdf_linear)
(
    double *restrict out,
    size_t length_y, size_t length_x_sort,
    double *restrict y,
    double *restrict x_sort,
);


#endif // L2SHM_CDF_H_
