#ifndef FTABLE_H
#define FTABLE_H

#include <gsl/gsl_spline.h>

#define Fsize (2000)
extern double Ftable[3][Fsize];

typedef struct FastPMFDInterp{
    size_t size;
    gsl_interp * F;
    gsl_interp * DF;
    // gsl_interp * DDF;
    gsl_interp_accel * acc;
} FastPMFDInterp;

// double getF(int F_id, double y);
void fastpm_fd_interp_init(FastPMFDInterp * FDinterp);
double fastpm_do_fd_interp(FastPMFDInterp * FDinterp, int F_id, double y);
void fastpm_fd_interp_destroy(FastPMFDInterp * FDinterp);

extern FastPMFDInterp *FDinterp;

#endif