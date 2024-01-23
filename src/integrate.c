//
// Created by Edward Parkinson on 23/01/2024.
//

#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

//
// Perform numerical integration on a given function which takes in a double and
// a void* of parameters.
//
double integrate(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound, double rtol) {
  double result;
  double error;
  size_t neval;

  gsl_function F;
  F.function = integrand;
  F.params = params;

  gsl_integration_workspace *qags_w = gsl_integration_workspace_alloc(1000);
  const int qags_status = gsl_integration_qags(&F, lower_bound, upper_bound, 0, rtol, 1000, qags_w, &result, &error);
  if (qags_status != GSL_SUCCESS) {
    if (qags_status == GSL_EROUND) {
      gsl_integration_workspace_free(qags_w);
      gsl_integration_romberg_workspace *romb_w = gsl_integration_romberg_alloc(30);
      const int romb_status = gsl_integration_romberg(&F, lower_bound, upper_bound, 0, rtol, &result, &neval, romb_w);

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
