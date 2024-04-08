#include <stdint.h>
#include <cblas.h>

#include <stdio.h>

#include <Rmath.h>


#include "l2shm.h"
#include "l2shm_random_R.h"


void L2SHM(random_combn)(size_t *restrict, size_t, size_t);


void L2SHM(select_batch)
(
    double *restrict U_batch,
    size_t n, size_t p, size_t b,
    double *restrict U,
    size_t *restrict indices
)
{
    L2SHM(random_combn)(indices, b, p-1);

    for(size_t ind = 0; ind < b; ++ind)
        cblas_dcopy(n, &U[indices[ind]*n], 1, &U_batch[ind*n], 1);

    return;
}

void L2SHM(random_combn)(size_t *restrict out, size_t length, size_t max){

  for(size_t ind = 0; ind < length; ++ind){
    size_t n = max-ind, s = SIZE_MAX - (SIZE_MAX % (n+1) + 1)*(SIZE_MAX % (n+1) < n), t;
    do{
      t = ((size_t) ((0x1.p+32)*unif_rand()) << 32) + (size_t) ((0x1.p+32)*unif_rand());
    } while(t > s);
    out[ind] = t%(n+1);
  }

  for(size_t ind1 = 1; ind1 < length; ++ind1){
    size_t ind2;
    for(ind2 = 0; out[ind2] <= out[ind1] && ind2 < ind1; ++ind2)
      ++out[ind1];
    size_t temp = out[ind1];
    for(size_t ind3 = ind1; ind3 > ind2; --ind3)
      out[ind3] = out[ind3-1];
    out[ind2] = temp;
  }

  return;
}
