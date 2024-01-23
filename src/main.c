//
// Created by Edward Parkinson on 22/01/2024.
//

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#include "atomic.h"
#include "python.h"

//
// Main function of the program
//
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
        const double recomb_sp = alpha_sp(&phot_top[xconfig[j].bfd_jump[k]], temperature, 0);
        Log("j %d: k=%d %g\n", j, k, recomb_sp);
      }
    }

    end_time = clock();
    const double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Time taken for cell %d: %f seconds\n", i, cpu_time_used);
  }

  return EXIT_SUCCESS;
}
