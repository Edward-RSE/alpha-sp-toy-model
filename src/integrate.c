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
// This uses Romberg integration.
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

//
// Perform numerical integration on a given function which takes in a double
// and a void * of parameters.
//
// This uses the new CQUAD integrator.
//
double integrate_cquad(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                       double rel_tol) {
  size_t n_evals = 0;
  double result = 0.0;
  double error = 0.0;

  gsl_function F;
  F.function = integrand;
  F.params = params;

  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(1000);
  gsl_integration_cquad(&F, lower_bound, upper_bound, 0, rel_tol, w, &result, &error, &n_evals);
  gsl_integration_cquad_workspace_free(w);

  return result;
}

//
// Perform numerical integration on a given function which takes in a double
// and a void * of parameters.
//
// This uses the QAG integrator.
//
double integrate_qag(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                     double rel_tol) {
  double result = 0.0;
  double error = 0.0;

  gsl_function F;
  F.function = integrand;
  F.params = params;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_integration_qag(&F, lower_bound, upper_bound, 0, rel_tol, 1000, GSL_INTEG_GAUSS31, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}


//
// Perform numerical integration on a given function which takes in a double
// and a void * of parameters.
//
// This uses a smaller QAGs integrator compared to the default.
//
double integrate_qags_small(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                            double rel_tol) {
  double result = 0.0;
  double error = 0.0;

  gsl_function F;
  F.function = integrand;
  F.params = params;

  const size_t limit = 100;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);
  gsl_integration_qags(&F, lower_bound, upper_bound, 0, rel_tol, limit, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

//
// Perform numerical integration on a given function which takes in a double
// and a void * of parameters.
//
// This performs the trapezium rule.
//
double integrate_trap(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                      double rel_tol) {
  const int max_iterations = 100;
  double result = 0.0;
  double prev_result = 0.0;

  for (int n = 1; n <= max_iterations; ++n) {
    double h = (upper_bound - lower_bound) / n;
    double sum = 0.5 * (integrand(lower_bound, params) + integrand(upper_bound, params));

    for (int i = 1; i < n; ++i) {
      double x = lower_bound + i * h;
      sum += integrand(x, params);
    }

    result = h * sum;

    if (n > 1 && fabs((result - prev_result) / result) < rel_tol) { break; }

    prev_result = result;
  }

  return result;
}
