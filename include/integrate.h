//
// Created by Edward Parkinson on 23/01/2024.
//

#ifndef NUM_INT_INTEGRATE_H
#define NUM_INT_INTEGRATE_H

double integrate(double (*integrand)(double, void *), void *params, double lower_bound, double upper_bound, double rtol);

#endif//NUM_INT_INTEGRATE_H
