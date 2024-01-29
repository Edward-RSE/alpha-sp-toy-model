//
// Created by Edward Parkinson on 23/01/2024.
//

#ifndef NUM_INT_INTEGRATE_H
#define NUM_INT_INTEGRATE_H

double integrate_default(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                         double rel_tol);
double integrate_romberg(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                         double rel_tol);
double integrate_cquad(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                       double rel_tol);
double integrate_qag(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                     double rel_tol);
double integrate_qags_small(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                            double rel_tol);
double integrate_trap(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                      double rel_tol);
double integrate_simp(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound,
                      double rel_tol);

#endif//NUM_INT_INTEGRATE_H
