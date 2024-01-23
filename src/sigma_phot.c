#include <stdio.h>
#include <math.h>

#include "atomic.h"
#include "python.h"



/* sigma_phot is used to calculate the photoionization cross-section */
double sigma_phot (struct topbase_phot *x_ptr, double freq)
{
  int nmax;
  double xsection;
  double frac;
  int nlast;

  if (freq < x_ptr->freq[0])
    return (0.0);               // Since this was below threshold
  if (freq == x_ptr->f)
    return (x_ptr->sigma);      // Avoid recalculating xsection

  if (x_ptr->nlast > -1)
  {
    nlast = x_ptr->nlast;
    if (x_ptr->freq[nlast] < freq && freq < x_ptr->freq[nlast + 1])
    {
      frac = (log (freq) - x_ptr->log_freq[nlast]) / (x_ptr->log_freq[nlast + 1] - x_ptr->log_freq[nlast]);
      xsection = exp ((1. - frac) * x_ptr->log_x[nlast] + frac * x_ptr->log_x[nlast + 1]);
      x_ptr->sigma = xsection;
      x_ptr->f = freq;

      return (xsection);
    }
  }

  /* Calculate the x-section */
  nmax = x_ptr->np;
  x_ptr->nlast = linterp (freq, &x_ptr->freq[0], &x_ptr->x[0], nmax, &xsection, 1);     // call linterp in log space

  x_ptr->sigma = xsection;
  x_ptr->f = freq;

  return (xsection);
}