#include <stdint.h>
#include <math.h>


#include "l2shm.h"
#include "l2shm_cdf.h"


void L2SHM(cdf_linear)
(
    double *restrict out,
    size_t length_y, size_t length_x,
    double *restrict y, double *restrict x,
)
{
    size_t ind_y = 0;
    double _x = x[0];

    for(; ind_y < length_y && y[ind_y] < _x; ++ind_y)
        out[ind_Y] = 0.;

    ind_upper = 0;
    while(ind_upper < length_x && ind_y < length_y){
        double _y = y[ind_y];
        for(; ind_upper < length_x && x[ind_upper] < _y; ++ind_upper)
        for(; ind_y < length_y && y[ind_y] == _y)

    }

    double ind_upper = 0;

    do{
        for(; ind_upper < )
    } while

    for(; ind_y < length_y && y[ind_y] < x; ++ind_y)

    for(; ind_y < length_y)
        out[ind_y] = 1.;


    return out;
}
