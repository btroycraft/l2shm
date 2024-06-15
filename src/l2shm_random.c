#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <stdint.h>
#include <cblas.h>
#include <math.h>


#include "l2shm.h"
#include "l2shm_random.h"


static size_t L2SHM(random_integer_weighted)
(
    size_t k,
    double *restrict Alpha
);

void L2SHM(random_batch)
(
    double *restrict U_batch,
    size_t n, size_t p, size_t b,
    double *restrict U,
    size_t *restrict indices
)
{
    L2SHM(random_combination)(indices, b, p-1);

    for(size_t ind = 0; ind < b; ++ind)
        cblas_dcopy(n, &U[indices[ind]*n], 1, &U_batch[ind*n], 1);

    return;
}

void L2SHM(random_uniform_sphere)
(
    double *restrict out,
    size_t dim
)
{
    GetRNGstate();
    for(size_t ind = 0; ind < dim; ++ind)
        out[ind] = norm_rand();
    PutRNGstate();

    cblas_dscal(dim, 1./cblas_dnrm2(dim, out, 1), out, 1);

    return;
}

void L2SHM(random_heat_sphere)
(
    double *restrict out,
    size_t dim,
    double *restrict mu, double T,
    size_t res
)
{
    double z = 1.;
    {
        double shape = (dim-2)/2.;

        double mult, scale;
        {
            double temp = 2.*T/res;
            mult = sqrt(temp);
            scale = 2.*temp;
        }

        GetRNGstate();
        for(size_t ind = 0; ind < res; ++ind){
            double temp = mult*norm_rand();
            z = (z+temp*sqrt(1.-z*z))/sqrt(1.+temp*temp+rgamma(shape, scale));
        }
    }

    for(size_t ind = 0; ind < dim; ++ind)
        out[ind] = norm_rand();
    PutRNGstate();

    {
        double nrm2 = cblas_dnrm2(dim, out, 1);
        double dot = cblas_ddot(dim, out, 1, mu, 1);

        double temp = sqrt((1.-z*z)/(nrm2*nrm2-dot*dot));
        cblas_dscal(dim, temp, out, 1);
        cblas_daxpy(dim, z-temp*dot, mu, 1, out, 1);
    }

    return;
}

void L2SHM(random_heat_sphere_mixture)
(
    double *restrict out,
    size_t dim, size_t k,
    double *restrict mu, double *restrict T, double *restrict Alpha,
    size_t res
)
{
    size_t index = L2SHM(random_integer_weighted)(k, Alpha);

    L2SHM(random_heat_sphere)(out, dim, &mu[index*dim], T[index], res);

    return;
}

void L2SHM(random_combination)
(
    size_t *restrict out,
    size_t length, size_t max
)
{
    GetRNGstate();
      for(size_t ind = 0; ind < length; ++ind){
        size_t n = max-ind, s = SIZE_MAX - (SIZE_MAX % (n+1) + 1)*(SIZE_MAX % (n+1) < n), t;
        do{
          t = ((size_t) ((0x1.p+32)*unif_rand()) << 32) + (size_t) ((0x1.p+32)*unif_rand());
        } while(t > s);
        out[ind] = t%(n+1);
      }
    PutRNGstate();

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

size_t L2SHM(random_integer)
(
    size_t max
)
{
    size_t s = SIZE_MAX - (SIZE_MAX % (max+1) + 1)*(SIZE_MAX % (max+1) < max), t;
    GetRNGstate();
    do{
      t = ((size_t) ((0x1.p+32)*unif_rand()) << 32) + (size_t) ((0x1.p+32)*unif_rand());
    } while(t > s);
    PutRNGstate();

    return t%(max+1);
}

static size_t L2SHM(random_integer_weighted)
(
    size_t k,
    double *restrict Alpha
)
{
    size_t out;
    {
        GetRNGstate();
        double val = unif_rand();
        PutRNGstate();
        for(out = 0; out < k-1; ++out){
            val -= Alpha[out];
            if(val < 0)
                break;
        }
    }

    return out;
}
