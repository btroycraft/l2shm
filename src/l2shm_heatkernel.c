#include <stdlib.h>
#include <math.h>

#include "l2shm.h"
#include "l2shm_heatkernel.h"


static inline double L2SHM(power_int)(double, size_t);


double L2SHM(heatkernel)(double x, double t, size_t dim, double terms)
{
    double d = dim;
    double dm4 = d-4., dm2i2 = 2./(d-2.);

    double t2 = t*t;

    double out = 1.;

    double p1 = 0., p2 = 1.;
    double ct = 1., cm = 1.;
    double mult = L2SHM(power_int)(t, dim-1);

    for(double ind = 1.; ind < terms; ++ind)
    {
        ct *= mult;
        cm += dm2i2;
        p1 += (dm4/ind+2.)*(p2*x-p1);

        out += ct*cm*p1;

        mult *= t2;
        ++ind;

        if(ind >= terms)
            break;

        ct *= mult;
        cm += dm2i2;
        p2 += (dm4/ind+2.)*(p1*x-p2);

        out += ct*cm*p2;

        mult *= t2;
    }

    return out;
}

double L2SHM(heatkernel_nox)(double t, size_t dim, double terms)
{
    double d = dim;
    double dm4 = d-4., dm2i2 = 2./(d-2.);

    double t2 = t*t;

    double out = 1.;

    double p1 = 0., p2 = 1.;
    double ct = 1., cm = 1.;
    double mult = L2SHM(power_int)(t, dim-1);

    for(double ind = 1.; ind < terms; ++ind)
    {
        ct *= mult;
        cm += dm2i2;
        p1 += (dm4/ind+2.)*(p2-p1);

        out += ct*cm*p1;

        mult *= t2;
        ++ind;

        if(ind >= terms)
            break;

        ct *= mult;
        cm += dm2i2;
        p2 += (dm4/ind+2.)*(p1-p2);

        out += ct*cm*p2;

        mult *= t2;
    }

    return out;
}

double L2SHM(heatkernel_dx)(double x, double t, size_t dim, double terms)
{
    double dm2 = dim-2.;

    double t2 = t*t;

    double out = 0.;

    double p1 = 0., p2 = 1.;
    double ct = 1., cm = dm2;
    double mult = L2SHM(power_int)(t, dim-1);

    for(double ind = 1.; ind < terms; ++ind)
    {
        ct *= mult;
        cm += 2.;
        p1 += (dm2/ind+2.)*(p2*x-p1);

        out += ct*cm*p2;

        mult *= t2;
        ++ind;

        if(ind >= terms)
            break;

        ct *= mult;
        cm += 2.;
        p2 += (dm2/ind+2.)*(p1*x-p2);

        out += ct*cm*p1;

        mult *= t2;
    }

    return out;
}

double L2SHM(heatkernel_dxx)(double x, double t, size_t dim, double terms)
{
    double dt2 = 2.*dim;
    double dm2 = dim-2.;

    double t2 = t*t;

    double out = 0.;

    double p1 = 0., p2 = 1.;
    double ct = 1., cm = dim*dm2;
    double mult = L2SHM(power_int)(t, dim-1);

    for(double ind = 1; ind < terms; ++ind)
    {
        ct *= mult;
        cm += dt2;

        out += ct*cm*p1;

        p1 += (dim/ind+2.)*(p2*x-p1);

        mult *= t2;
        ++ind;

        if(ind >= terms)
            break;

        ct *= mult;
        cm += dt2;

        out += ct*cm*p2;

        p2 += (dim/ind+2.)*(p1*x-p2);

        mult *= t2;
    }

    return out;
}

double L2SHM(heatkernel_dt)(double x, double t, size_t dim, double terms)
{
    double dm4 = dim-4.;
    double dm2i12 = 12./(dim-2.);
    double t2 = t*t;

    double out = 0.;

    double p1 = 0., p2 = 1.;
    double ct = 1./t, cm = 0.;
    double add1 = 6., add2 = 2./(dim-2.)+dim+1.;
    double mult = L2SHM(power_int)(t, dim-1);

    for(double ind = 1.; ind < terms; ++ind)
    {
        ct *= mult;
        cm += add2;
        p1 += (dm4/ind+2.)*(p2*x-p1);

        out += ct*cm*p1;

        add1 += dm2i12;
        add2 += add1;
        mult *= t2;

        ++ind;

        if(ind >= terms)
            break;

        ct *= mult;
        cm += add2;
        p2 += (dm4/ind+2.)*(p1*x-p2);

        out += ct*cm*p2;

        add1 += dm2i12;
        add2 += add1;
        mult *= t2;
    }

    return out;
}

void L2SHM(heatkernel_combn_dxdt)
(
    double * out, double * out_dx, double * out_dt,
    double x, double t, size_t dim, double terms
)
{
    double d = dim;
    double dm2 = d-2., dm4 = d-4.;
    double dm2i2 = 2./dm2;

    double t2 = t*t, ti = 1./t;
    double ti2 = 2.*ti;

    double _out = 1., _out_dx = 0., _out_dt = 0.;

    double p11 = 0., p12 = 1., p21 = 0., p22 = 1.;
    double ct = 1., cm = 1., cd = 0.;
    double add = (d-1.)*ti;

    double mult = L2SHM(power_int)(t, dim-1);

    for(double ind = 1.; ind < terms; ++ind)
    {
        ct *= mult;
        cm += dm2i2;
        cd += add;

        {
            double temp = 1./ind;
            p11 += (dm4*temp+2.)*(p12*x-p11);
            p21 += (dm2*temp+2.)*(p22*x-p21);
        }
        {
            double temp1 = ct*cm;
            _out_dx += temp1*p22;
            double temp2 = temp1*p11;
            _out += temp2;
            _out_dt += cd*temp2;
        }

        add += ti2;
        mult *= t2;
        ++ind;

        if(ind >= terms)
            break;

        ct *= mult;
        cm += dm2i2;
        cd += add;

        {
            double temp = 1./ind;
            p12 += (dm4*temp+2.)*(p11*x-p12);
            p22 += (dm2*temp+2.)*(p21*x-p22);
        }
        {
            double temp1 = ct*cm;
            _out_dx += temp1*p21;
            double temp2 = temp1*p12;
            _out += temp2;
            _out_dt += cd*temp2;
        }

        add += ti2;
        mult *= t2;
    }
    _out_dx *= dm2;

    *out = _out;
    *out_dx = _out_dx;
    *out_dt = _out_dt;

    return;
}

void L2SHM(heatkernel_combn_noxdt)
(
    double * out, double * out_dt,
    double t, size_t dim, double terms
)
{
    double d = dim;
    double dm4 = d-4., dm2i2 = 2./(d-2.);

    double t2 = t*t, ti = 1./t;
    double ti2 = 2.*ti;

    double _out = 1., _out_dt = 0.;

    double p1 = 0., p2 = 1.;
    double ct = 1., cm = 1., cd = 0.;
    double add = (d-1.)*ti;

    double mult = L2SHM(power_int)(t, dim-1);

    for(double ind = 1.; ind < terms; ++ind)
    {
        ct *= mult;
        cm += dm2i2;
        cd += add;

        p1 += (dm4/ind+2.)*(p2-p1);

        {
            double temp = ct*cm*p1;
            _out += temp;
            _out_dt += cd*temp;
        }

        add += ti2;
        mult *= t2;
        ++ind;

        if(ind >= terms)
            break;

        ct *= mult;
        cm += dm2i2;
        cd += add;

        p2 += (dm4/ind+2.)*(p1-p2);

        {
            double temp = ct*cm*p2;
            _out += temp;
            _out_dt += cd*temp;
        }

        add += ti2;
        mult *= t2;
    }

    *out = _out;
    *out_dt = _out_dt;

    return;
}

double L2SHM(heatkernel_test)(double x, double t, size_t dim, double terms)
{
    double dm2 = dim-2.;
    double dm4 = dim-4.;
    double dm2i2 = 2./(dim-2.);

    double out = 1.;

    double p1 = 0., p2 = 1.;

    for(double ind = 1.; ind < terms; ++ind)
    {
        p1 += (dm4/ind+2.)*(p2*x-p1);

        out += pow(t, ind*(ind+dm2))*(dm2i2*ind+1.)*p1;

        ++ind;

        if(ind >= terms)
            break;

        p2 += (dm4/ind+2.)*(p1*x-p2);

        out += pow(t, ind*(ind+dm2))*(dm2i2*ind+1.)*p2;
    }

    return out;
}

double L2SHM(heatkernel_dx_test)(double x, double t, size_t dim, double terms)
{
    double dm2 = dim-2.;

    double out = 0.;

    double p1 = 0., p2 = 1.;

    for(double ind = 1; ind < terms; ++ind)
    {
        p1 += (dm2/ind+2.)*(p2*x-p1);

        out += pow(t, ind*(ind+dm2))*(2.*ind+dm2)*p2;

        ++ind;

        if(ind >= terms)
            break;

        p2 += (dm2/ind+2.)*(p1*x-p2);

        out += pow(t, ind*(ind+dm2))*(2.*ind+dm2)*p1;
    }

    return out;
}

double L2SHM(heatkernel_dxx_test)(double x, double t, size_t dim, double terms)
{
    double dt2 = 2.*dim;
    double dm2 = dim-2.;
    double dtdm2 = dim*dm2;

    double out = 0.;

    double p1 = 0., p2 = 1.;

    for(double ind = 1; ind < terms; ++ind)
    {
        out += pow(t, ind*(ind+dm2))*(dt2*ind+dtdm2)*p1;

        p1 += (dim/ind+2.)*(p2*x-p1);

        ++ind;

        if(ind >= terms)
            break;

        out += pow(t, ind*(ind+dm2))*(dt2*ind+dtdm2)*p2;

        p2 += (dim/ind+2.)*(p1*x-p2);
    }

    return out;
}

double L2SHM(heatkernel_dt_test)(double x, double t, size_t dim, double terms)
{
    double dm2 = dim-2.;
    double dm4 = dim-4.;
    double dm2i2 = 2./(dim-2.);

    double out = 0.;

    double p1 = 0., p2 = 1.;

    for(double ind = 1.; ind < terms; ++ind)
    {
        p1 += (dm4/ind+2.)*(p2*x-p1);

        out += pow(t, ind*(ind+dm2)-1.)*(dm2i2*ind+1.)*ind*(ind+dm2)*p1;

        ++ind;

        if(ind >= terms)
            break;

        p2 += (dm4/ind+2.)*(p1*x-p2);

        out += pow(t, ind*(ind+dm2)-1.)*(dm2i2*ind+1.)*ind*(ind+dm2)*p2;
    }

    return out;
}


static inline double L2SHM(power_int)(double a, size_t power)
{
    double out = 1.;
    double a2k = a;
    for(size_t bits = power; bits > 0; bits >>= 1)
    {
        if(bits & 1)
            out *= a2k;
        a2k *= a2k;
    }
    return out;
}
