/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  ngenic.h
 *
 *  \brief definition of a class for the construction of cosmological initial conditions
 */

#ifndef NGENIC_H
#define NGENIC_H

#ifdef NGENIC

#ifndef PERIODIC
#error NGENIC requires PERIODIC
#endif

#include <fftw3.h>

#ifdef DOUBLEPRECISION_FFTW
typedef double fft_real;
typedef fftw_complex fft_complex;
#else
typedef float fft_real;
typedef fftwf_complex fft_complex;
#endif

#include <gsl/gsl_odeiv2.h>
#include "../data/simparticles.h"
#include "../pm/pm_mpi_fft.h"
#include "../data/mymalloc.h"

class ngenic : public pm_mpi_fft
{
 private:
  simparticles *Sp;

 public:
  ngenic(MPI_Comm comm, simparticles *Sp_ptr) : setcomm(comm), pm_mpi_fft(comm) /* constructor */ { Sp = Sp_ptr; }

 public:
  void ngenic_displace_particles(void);

  void create_grid(void);

 private:
  double ngenic_power_spec(double k);
  double ngenic_f1_omega(double a);
  double ngenic_f2_omega(double a);
  double ngenic_Phi2fac(double a);
  double ngenic_growth_factor(double astart, double aend);
  void ngenic_initialize_powerspectrum(void);
  void free_power_table(void);

  double Dplus;

  unsigned int *seedtable;

  fft_plan myplan;
  size_t maxfftsize;

  struct partbuf
  {
    MyIntPosType IntPos[3];
  };
  partbuf *partin, *partout;

  size_t nimport, nexport;

  size_t *Sndpm_count, *Sndpm_offset;
  size_t *Rcvpm_count, *Rcvpm_offset;

  gsl_rng *rnd_generator;
  gsl_rng *rnd_generator_conjugate;

  struct disp_data
  {
    fft_real deltapos[3];
  };

  disp_data *Pdisp;

  void ngenic_distribute_particles();
  void ngenic_setup_modes_in_kspace(fft_complex *fft_of_grid);
  void ngenic_readout_disp(fft_real *grid, int axis, double pfac, double vfac);
  void ngenic_initialize_ffts(void);
  void ngenic_get_derivate_from_fourier_field(int axes1, int axes2, fft_complex *fft_of_grid);
  void ngenic_compute_transform_of_source_potential(fft_real *pot);
  void print_spec(void);

  double R8;

  double AA, BB, CC;
  double nu;
  double Norm;

  int NPowerTable;

  struct pow_table
  {
    double logk, logD;
    bool operator<(const pow_table &other) const { return logk < other.logk; }
  };
  pow_table *PowerTable;

  double ngenic_powerspec_tabulated(double k);
  double ngenic_powerspec_efstathiou(double k);
  double ngenic_powerspec_eh(double k);
  double ngenic_tophat_sigma2(double R);
  double ngenic_tk_eh(double k);
  double ngenic_growth(double a);
  void read_power_table(void);

  static double ngenic_growth_int(double a, void *param)
  {
    return pow(a / (All.Omega0 + (1 - All.Omega0 - All.OmegaLambda) * a + All.OmegaLambda * a * a * a), 1.5);
  }

  double fnl(double x, double A, double B, double alpha, double beta, double V, double gf) /* Peacock & Dodds formula */
  {
    return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) / (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)),
                   1 / beta);
  }
/////////////////////////// Chen's modification ////////////////////////  
/** Exact solution is modified accroding to the FastPM "libfastpm/cosmology.c" file. */
typedef struct {
    double y0;
    double y1;
    double y2;
    double y3;
} ode_soln;

static int growth_ode(double a, const double y[], double dyda[], void* params)
{
    
    const double E_ct = Driftfac.hubble_function(a)  / All.Hubble;
    const double dEda = Driftfac.DHubbleEaDa(a);
    double Omega_sourc_cdm = All.Omega0 / (a * a * a) / E_ct / E_ct;
    double dydlna[4];
    dydlna[0] = y[1];
    dydlna[1] = - (2. + a / E_ct * dEda) * y[1] + 1.5 * Omega_sourc_cdm * y[0];
    dydlna[2] = y[3];
    dydlna[3] = - (2. + a / E_ct * dEda) * y[3] + 1.5 * Omega_sourc_cdm * (y[2] - y[0]*y[0]);
    
    //divide by  a to get dyda
    for (int i=0; i<4; i++){
        dyda[i] = dydlna[i] / a;
    }
    
    return GSL_SUCCESS;
}

static ode_soln growth_ode_solve(double a)
{
    /* This returns an array of {d1, F1, d2, F2} (unnormalised) */
    gsl_odeiv2_system F;
    F.function = &growth_ode;
    F.jacobian = NULL;
    F.dimension = 4;

    gsl_odeiv2_driver * drive 
        = gsl_odeiv2_driver_alloc_standard_new(&F,
                                               gsl_odeiv2_step_rkf45, 
                                               1e-6,
                                               1e-8,
                                               1e-8,
                                               1,
                                               1);
    
    // assume matter domination
    double aini = 0.00625;  // FIXME: need to make sure this is less than the starting a. For now using z=159.
    double yini[4];
    yini[0] = aini;
    yini[1] = aini;
    yini[2] = - 3./7. * aini*aini;
    yini[3] = 2 * yini[2];
    
    int status = gsl_odeiv2_driver_apply(drive, &aini, a, yini);
    if (status != GSL_SUCCESS) {
        printf("Growth ODE unsuccesful at a=%g.", a);
    }
    
    gsl_odeiv2_driver_free(drive);
    
    ode_soln soln;
    soln.y0 = yini[0];
    soln.y1 = yini[1];
    soln.y2 = yini[2];
    soln.y3 = yini[3];
    
    return soln;
}
/////////////////////////// Chen's modification ////////////////////////  

  struct myparams
  {
    double R;
    ngenic *obj;
  };

  static double sigma2_int(double lnk, void *param)
  {
    myparams *par   = (myparams *)param;
    double r_tophat = par->R;

    double k   = exp(lnk);
    double kr  = r_tophat * k;
    double kr2 = kr * kr;
    double kr3 = kr2 * kr;

    if(kr < 1e-8)
      return 0;

    double w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
    double x = 4 * M_PI * k * k * w * w * par->obj->ngenic_power_spec(k);

    return k * x;
  }
};

#endif
#endif
