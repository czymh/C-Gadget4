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
#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/dtypes.h"
#include "../main/main.h"
#include "gadgetconfig.h"

/*Stefan-Boltzmann constant in cgs units*/
#define  STEFAN_BOLTZMANN 5.670373e-5

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
    return ( 3*All.DE_wa - (1 + All.DE_w0 + All.DE_wa)/a ) * OmegaLTimesHubbleSquare(a);
  }

  static double Omega_ncdmTimesHubbleSquare(double a)
  {
    return All.M_nu_all / 93.14 / (All.HubbleParam*All.HubbleParam) / (a * a * a);
  }

  static double Omega_RadiationTimesHubbleSquare(double a)
  {
    /* Consider the Radiation: Omega_g = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H^2) */
    double Omega_Gamma = 4 * STEFAN_BOLTZMANN
                           * pow(All.CMBTemperature, 4)
                           * (8 * M_PI * GRAVITY)
                           / (3*pow(CLIGHT, 3)*HUBBLE*HUBBLE)
                           / (All.HubbleParam*All.HubbleParam);
    /* the ultra-relativistiv parts are treated as radiation*/
    double fur = 7/8*pow((4/11), (4/3))*All.Nur;
    return (1+fur) * Omega_Gamma / (a * a * a * a);
  }

  /** Accually, hubble_function = E(a) = H(a)/h0 */
  static double hubble_function(double a)
  {
    /** Omega0 = Omega_CDM + Omega_Baryon*/
    double hubble_a = All.Omega0 / (a * a * a);
    /* The Dynamic Dark Energy, equation of state:  w(a) = w0 + wa(1 -a)  */
    hubble_a += OmegaLTimesHubbleSquare(a);
    /* Neutrino is only considered as the smooth component in the background  */
    /* only for the degenerate neutrino mass */
    hubble_a += Omega_ncdmTimesHubbleSquare(a);
    double Omeganu = Omega_ncdmTimesHubbleSquare(1.0);
    /* All radiations including ultra-relativistic components (e.g. massless neutrinos)*/
    hubble_a += Omega_RadiationTimesHubbleSquare(a);
    /** modified curvature */
    hubble_a += (1 - All.Omega0 - All.OmegaLambda - Omeganu) / (a * a);
    hubble_a = All.Hubble * sqrt(hubble_a);
    return hubble_a;
  }

  double DHubbleEaDa(double a)
  {
    double H_ct = hubble_function(a) / All.Hubble;
    double Omeganu = Omega_ncdmTimesHubbleSquare(1.0);
    return 0.5 / H_ct * (- 4 * Omega_RadiationTimesHubbleSquare(a) / a
                         - 3 * All.Omega0 / (a * a * a * a)
                         - 2 * (1 - All.Omega0 - All.OmegaLambda - Omeganu) / (a * a * a)
                         - 4 * Omega_ncdmTimesHubbleSquare(a) / a
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
