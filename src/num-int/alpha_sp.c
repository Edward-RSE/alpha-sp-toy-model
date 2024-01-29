//
// Created by Edward Parkinson on 23/01/2024.
//

#include <math.h>
#include <stdio.h>

#include "atomic.h"
#include "integrate.h"
#include "python.h"

//
// `sigma_phot` is used to calculate the photoionization cross-section for a
// given topbase photoionization entry and frequency. This function is taken
// from Python
//
double sigma_phot(struct topbase_phot *x_ptr, double freq) {
  int nmax;
  double xsection;
  double frac;
  int nlast;

  if (freq < x_ptr->freq[0]) {
    return (0.0);// Since this was below threshold
  }

  if (freq == x_ptr->f) {
    return (x_ptr->sigma);// Avoid recalculating xsection
  }

  if (x_ptr->nlast > -1) {
    nlast = x_ptr->nlast;
    if (x_ptr->freq[nlast] < freq && freq < x_ptr->freq[nlast + 1]) {
      frac = (log(freq) - x_ptr->log_freq[nlast]) / (x_ptr->log_freq[nlast + 1] - x_ptr->log_freq[nlast]);
      xsection = exp((1. - frac) * x_ptr->log_x[nlast] + frac * x_ptr->log_x[nlast + 1]);
      x_ptr->sigma = xsection;
      x_ptr->f = freq;

      return (xsection);
    }
  }

  /* Calculate the x-section */
  nmax = x_ptr->np;
  x_ptr->nlast = linterp(freq, &x_ptr->freq[0], &x_ptr->x[0], nmax, &xsection, 1);// call linterp in log space
  x_ptr->sigma = xsection;
  x_ptr->f = freq;

  return (xsection);
}

//
// This is the struct we'll use to pass the integration parameters to the GSL
// numerical routines
//
struct integration_parameters {
  double temperature;
  double freq_lower;
  struct topbase_phot *phot;
};

//
// Returns the value of the integrand in the calculation for the spontaneous
// recombination coefficient
//
double alpha_sp_integration(double freq, void *params) {
  const struct integration_parameters *p = (struct integration_parameters *) params;
  const double temperature = p->temperature;
  const double freq_lower = p->freq_lower;
  struct topbase_phot *phot = p->phot;

  if (freq < freq_lower) { return 0.0; }

  const double x_section = sigma_phot(phot, freq);
  double integrand = x_section * freq * freq * exp(H_OVER_K * (freq_lower - freq) / temperature);

  return integrand;
}

//
// Calculate the spontaneous recombination coefficient for a given temperature
// and photoionization level
//
#define ALPHA_SP_CONSTANT 5.79618e-36
double alpha_sp(struct topbase_phot *phot, const double temperature, int mode,
                double (*integrator)(double (*integrand)(double, void *), void *, double, double, double)) {
  const double rtol = 1e-4;// hardcoded to 1e-4, like in Python
  const double freq_lower = phot->freq[0];
  double freq_upper = phot->freq[phot->np - 1];

  if ((H_OVER_K * (freq_upper - freq_lower) / temperature) > ALPHA_MATOM_NUMAX_LIMIT) {
    freq_upper = freq_lower + temperature * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }

  struct integration_parameters params = {temperature, freq_lower, phot};
  double recomb_sp_value = integrator(alpha_sp_integration, &params, freq_lower, freq_upper, rtol);

  if (phot->macro_info == TRUE && geo.macro_simple == FALSE) {
    recomb_sp_value *= xconfig[phot->nlev].g / xconfig[phot->uplev].g * pow(temperature, -1.5);
  } else {
    recomb_sp_value *= xconfig[phot->nlev].g / xconfig[phot->nion + 1].g * pow(temperature, -1.5);
  }

  recomb_sp_value *= ALPHA_SP_CONSTANT;

  return recomb_sp_value;
}
