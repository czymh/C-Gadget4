/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  driftfac.h
 *
 *  \brief declares a class for supporting cosmological drift/kick factor calculations
 */

#ifndef DRIFTFAC_H
#define DRIFTFAC_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/dtypes.h"
#include "../main/main.h"
#include "gadgetconfig.h"
#include "./Ftable.h"

/*Stefan-Boltzmann constant in cgs units*/
#define STEFAN_BOLTZMANN 5.670373e-5
#define kB 8.617333262145e-5   // boltzman in eV/K

class driftfac
{
 public:
  void init_drift_table(void);
  double get_drift_factor(integertime time0, integertime time1);
  double get_gravkick_factor(integertime time0, integertime time1);
  double get_hydrokick_factor(integertime time0, integertime time1);
  double get_comoving_distance(integertime time0);
  double get_comoving_distance_for_scalefactor(double ascale);
  double get_scalefactor_for_comoving_distance(double dist);
  integertime get_gravkick_factor_inverse(double fac);

  static double OmegaLTimesHubbleSquare(double a)
  {
    double exponent = (a-1)*All.DE_wa-(1 + All.DE_w0 + All.DE_wa)*log(a);
    return All.OmegaLambda * exp(3 * exponent); 
  }
  static double DOmegaLTimesHubbleSquareDa(double a)
  {
    return ( 3*All.DE_wa - 3*(1 + All.DE_w0 + All.DE_wa)/a ) * OmegaLTimesHubbleSquare(a);
  }
  
  static double Gamma_nu(void)
  {
    /* nu to photon temp ratio today */
    /** Added Parameter Nncdm*/
    /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    // This is consistent with CLASS but different to FastPM defs.
    // double N_eff = All.Nur + All.Nncdm*1.0132 
    return pow(4./11., 1./3.);
  }

  static double OmegaGamma(void)
  {
    /* Consider the Radiation: Omega_g = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H^2) */
    return 4 * STEFAN_BOLTZMANN
             * pow(All.CMBTemperature, 4)
             * (8 * M_PI * GRAVITY)
             / (3*pow(CLIGHT, 3)*HUBBLE*HUBBLE)
             / (All.HubbleParam*All.HubbleParam);
  } 

  static double Fconst(int ncdm_id)
  {
      /* This is a cosmology dependent constant 
         which is the argument divided by a of F, DF, DDF */
      double T_nu = Gamma_nu()*All.CMBTemperature;
      double m_nu_list[3] = {All.m_nu1, All.m_nu2, All.m_nu3};
      return m_nu_list[ncdm_id] / (kB * T_nu);
  }

  static double getFtable(int F_id, double y)
  {
    /* Gets the interpolated value of Ftable[F_id] at y
       F_id: 1 for F, 2 for F', 3 for F'' */
    if (y > 0) {
        // return getF(F_id, y);
        return fastpm_do_fd_interp(FDinterp, F_id, y);
    }else{
        return 0; 
    }
  }
 
  static double Omega_ncdmTimesHubbleSquare(double a)
  {
    /* Omega_ncdm(a) * E(a)^2 */
    // return All.M_nu_all / 93.14 / (All.HubbleParam*All.HubbleParam) / (a * a * a); // only for the CDM-like
    double A = 15. / pow(M_PI, 4) * pow(Gamma_nu(), 4) * OmegaGamma();     
    double Fall = 0;
    for (int i=0; i<All.Nncdm; i++) {
        double Fci = Fconst(i); 
        Fall += getFtable(1, Fci*a); //row 1 for F
    }
    return A / (a*a*a*a) * Fall;
  }
  
  static double DOmega_ncdmTimesHubbleSquareDa(double a)
  {
    double A = 15. / pow(M_PI, 4) * pow(Gamma_nu(), 4) * OmegaGamma();
    double OncdmESq = Omega_ncdmTimesHubbleSquare(a);
    
    double FcDF = 0;
    for (int i=0; i<All.Nncdm; i++) {
        double Fc = Fconst(i);
        double DF = getFtable(2, Fc*a);
        FcDF += Fc * DF;    //row 2 for F'
    }
    
    return -4. / a * OncdmESq + A / (a*a*a*a) * FcDF;
  }

  static double Omega_RadiationTimesHubbleSquare(double a)
  {
    double Omega_Gamma = OmegaGamma();
    /* the ultra-relativistiv parts are treated as radiation*/
    double fur = 7./8. * pow(Gamma_nu(), 4) * All.Nur;
    // double fur = 7./8. * pow(4./11., 1./3.) * All.Nur
    return (1+fur) * Omega_Gamma / (a * a * a * a);
  }

  /** Accually, hubble_function = E(a) = H(a)/h0 */
  static double hubble_function(double a)
  {
    // if (All.Nncdm > 0){
    //   fastpm_fd_interp_init(FDinterp);
    // } /** initialize the Ftable integration */
    /** Omega0 = Omega_CDM + Omega_Baryon*/
    double hubble_a = All.Omega0 / (a * a * a);
    /* The Dynamic Dark Energy, equation of state:  w(a) = w0 + wa(1 -a)  */
    hubble_a += OmegaLTimesHubbleSquare(a);
    /* Neutrino is only considered as the smooth component in the background  */
    /* only for the degenerate neutrino mass */
    hubble_a += Omega_ncdmTimesHubbleSquare(a);
    //double Omeganu = Omega_ncdmTimesHubbleSquare(1.0);
    /* All radiations including ultra-relativistic components (e.g. massless neutrinos)*/
    hubble_a += Omega_RadiationTimesHubbleSquare(a);
    /** modified curvature */
    //double Omegak = 1 - All.Omega0 - All.OmegaLambda - Omeganu - Omega_RadiationTimesHubbleSquare(1.0)
    double Omegak = 0; // Set 0 
    hubble_a += Omegak / (a * a);
    hubble_a = All.Hubble * sqrt(hubble_a);
    // if (FDinterp){
    //   fastpm_fd_interp_init(FDinterp);
    // } /** free memory for stable*/
    return hubble_a;
  }

  double DHubbleEaDa(double a)
  {
    double H_ct = hubble_function(a) / All.Hubble;
    //double Omeganu = Omega_ncdmTimesHubbleSquare(1.0);
    double Omegak = 0; // 1 - All.Omega0 - All.OmegaLambda - Omeganu - Omega_RadiationTimesHubbleSquare(1.0)
    return 0.5 / H_ct * (- 4 * Omega_RadiationTimesHubbleSquare(a) / a
                         - 3 * All.Omega0 / (a * a * a * a)
                         - 2 * Omegak / (a * a * a)
                         + DOmega_ncdmTimesHubbleSquareDa(a)
                         + DOmegaLTimesHubbleSquareDa(a)
                        );
  }

 private:
#define DRIFT_TABLE_LENGTH 1000

  /** table for the cosmological drift factors */
  double DriftTable[DRIFT_TABLE_LENGTH];

  /** table for the cosmological kick factor for gravitational forces */
  double GravKickTable[DRIFT_TABLE_LENGTH];

  /** table for the cosmological kick factor for hydrodynmical forces */
  double HydroKickTable[DRIFT_TABLE_LENGTH];

  double logTimeBegin;
  double logTimeMax;

  static double drift_integ(double a, void *param)
  {
    double h = hubble_function(a);

    return 1 / (h * a * a * a);
  }

  static double gravkick_integ(double a, void *param)
  {
    double h = hubble_function(a);

    return 1 / (h * a * a);
  }

  static double hydrokick_integ(double a, void *param)
  {
    double h = hubble_function(a);

    return 1 / (h * pow(a, 3 * GAMMA_MINUS1) * a);
  }
};

extern driftfac Driftfac;

#endif
