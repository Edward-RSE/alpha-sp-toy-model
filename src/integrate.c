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
// This implementation is what it more or less is in Python, minus some changes
// to how integration parameters are passed.
//
double integrate_default(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                         double rel_tol) {
  double result;
  double error;
  size_t neval;

  gsl_function F;
  F.function = integrand;
  F.params = params;

  gsl_integration_workspace *qags_w = gsl_integration_workspace_alloc(1000);
  const int qags_status = gsl_integration_qags(&F, lower_bound, upper_bound, 0, rel_tol, 1000, qags_w, &result, &error);
  if (qags_status != GSL_SUCCESS) {
    // fallback to romberg integration when something goes wrong
    if (qags_status == GSL_EROUND) {
      gsl_integration_workspace_free(qags_w);
      gsl_integration_romberg_workspace *romb_w = gsl_integration_romberg_alloc(30);
      const int romb_status =
          gsl_integration_romberg(&F, lower_bound, upper_bound, 0, rel_tol, &result, &neval, romb_w);

      if (romb_status != GSL_SUCCESS) { fprintf(stderr, "numerical integration error\n"); }
      gsl_integration_romberg_free(romb_w);
    }
  } else {
    gsl_integration_workspace_free(qags_w);
  }

  return result;
}

//
// Perform numerical integration on a given function which takes in a double and
// a void* of parameters.
//
// This is a new function we will use to experiment with.
//
double integrate_romberg(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                     double rel_tol) {
  size_t n_evals = 0;
  double result = 0.0;

  gsl_function F;
  F.function = integrand;
  F.params = params;

  gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc(30);
  gsl_integration_romberg(&F, lower_bound, upper_bound, 0, rel_tol, &result, &n_evals, w);
  gsl_integration_romberg_free(w);

  return result;
}
