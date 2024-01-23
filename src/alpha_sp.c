//
// Created by Edward Parkinson on 23/01/2024.
//

#include <math.h>
#include <stdio.h>

#include "atomic.h"
#include "python.h"
#include "integrate.h"

#define ALPHA_SP_CONSTANT 5.79618e-36

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

  if (freq < freq_lower) {
    return 0.0;
  }

  const double x_section = sigma_phot(phot, freq);
  double integrand = x_section * freq * freq * exp(H_OVER_K * (freq_lower - freq) / temperature);

  return integrand;
}

//
// Calculate the spontaenous recombination coefficient for a given temperature
// and photoionization level
//
double alpha_sp(struct topbase_phot *phot, const double temperature, int mode) {
  const double rtol = 1e-4;// hardcoded to 1e-4, like in Python
  const double freq_lower = phot->freq[0];
  double freq_upper = phot->freq[phot->np - 1];

  if ((H_OVER_K * (freq_upper - freq_lower) / temperature) > ALPHA_MATOM_NUMAX_LIMIT) {
    freq_upper = freq_lower + temperature * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }

  struct integration_parameters params = {temperature, freq_lower, phot};
  double recomb_sp_value = integrate(alpha_sp_integration, &params, freq_lower, freq_upper, rtol);

  if (phot->macro_info == TRUE && geo.macro_simple == FALSE) {
    recomb_sp_value *= xconfig[phot->nlev].g / xconfig[phot->uplev].g * pow(temperature, -1.5);
  } else {
    recomb_sp_value *= xconfig[phot->nlev].g / xconfig[phot->nion + 1].g * pow(temperature, -1.5);
  }

  recomb_sp_value *= ALPHA_SP_CONSTANT;

  return recomb_sp_value;
}
