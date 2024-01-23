#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "atomic.h"
#include "python.h"

#define ALPHA_SP_CONSTANT 5.79618e-36

struct recomb_sp_parameters {
  double temperature;
  double freq_lower;
  struct topbase_phot *pi;
};

double recomb_sp_integrand(double freq, void *params) {
  const struct recomb_sp_parameters *p = (struct recomb_sp_parameters *) params;
  const double temperature = p->temperature;
  const double freq_lower = p->freq_lower;
  struct topbase_phot *phot = p->pi;

  if (freq < freq_lower) {
    return 0.0;
  }

  const double x_section = sigma_phot(phot, freq);
  double integrand = x_section * freq * freq * exp(H_OVER_K * (freq_lower - freq) / temperature);

  return integrand;
}

double integrate(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound) {
  double result;
  double error;
  size_t neval;

  gsl_function F;
  F.function = integrand;
  F.params = params;

  const double EPS = 1e-4;

  gsl_integration_workspace *qags_w = gsl_integration_workspace_alloc(1000);
  const int qags_status = gsl_integration_qags(&F, lower_bound, upper_bound, 0, EPS, 1000, qags_w, &result, &error);
  if (qags_status != GSL_SUCCESS) {
    if (qags_status == GSL_EROUND) {
      fprintf(stdout, "GSL_EROUND\n");
      gsl_integration_workspace_free(qags_w);
      gsl_integration_romberg_workspace *romb_w = gsl_integration_romberg_alloc(30);
      const int romb_status = gsl_integration_romberg(&F, lower_bound, upper_bound, 0, EPS, &result, &neval, romb_w);
      if (romb_status != GSL_SUCCESS) {
        fprintf(stderr, "numerical integration error\n");
      }
      gsl_integration_romberg_free(romb_w);
    } else {
      fprintf(stderr, "numerical integration error\n");
    }
  } else {
    gsl_integration_workspace_free(qags_w);
  }

  return result;
}

double get_recomb_sp(struct topbase_phot *phot, const double temperature) {
  const double freq_lower = phot->freq[0];
  double freq_upper = phot->freq[phot->np - 1];

  if ((H_OVER_K * (freq_upper - freq_lower) / temperature) > ALPHA_MATOM_NUMAX_LIMIT) {
    freq_upper = freq_lower + temperature * ALPHA_MATOM_NUMAX_LIMIT / H_OVER_K;
  }

  struct recomb_sp_parameters params = {temperature, freq_lower, phot};
  double recomb_sp_value = integrate(recomb_sp_integrand, &params, freq_lower, freq_upper);

  if (phot->macro_info == TRUE && geo.macro_simple == FALSE) {
    recomb_sp_value *= xconfig[phot->nlev].g / xconfig[phot->uplev].g * pow(temperature, -1.5);
  } else {
    recomb_sp_value *= xconfig[phot->nlev].g / xconfig[phot->nion + 1].g * pow(temperature, -1.5);
  }

  recomb_sp_value *= ALPHA_SP_CONSTANT;

  return recomb_sp_value;
}


int main(int argc, char **argv) {
  clock_t start_time, end_time;

  geo.ioniz_mode = 9;
  Log_set_verbosity(0);
  get_atomic_data("data/h10_hetop_standard80.dat");
  gsl_set_error_handler_off();

  const int num_cells = 1;
  const double temperature = 40000;

  for (int i = 0; i < num_cells; ++i) {
    start_time = clock();

    for (int j = 0; j < nlevels_macro; ++j) {
      for (int k = 0; k < xconfig[j].n_bfd_jump; ++k) {
        const double recomb_sp = get_recomb_sp(&phot_top[xconfig[j].bfd_jump[k]], temperature);
        Log("j %d: k=%d %g\n", j, k, recomb_sp);
      }
    }

    end_time = clock();
    const double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Time taken for cell %d: %f seconds\n", i, cpu_time_used);
  }

  return EXIT_SUCCESS;
}
