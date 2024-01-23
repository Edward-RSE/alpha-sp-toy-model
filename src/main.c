//
// Created by Edward Parkinson on 22/01/2024.
//

#include <gsl/gsl_errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "atomic.h"
#include "integrate.h"
#include "python.h"

//
// Print a divider for model initialisation
//
void print_initialise_divider(void) {
  printf("---------------------------------------\n");
  printf(" ___       _ _   _       _ _          \n");
  printf("|_ _|_ __ (_) |_(_) __ _| (_)___  ___ \n");
  printf(" | || '_ \\| | __| |/ _` | | / __|/ _ \\\n");
  printf(" | || | | | | |_| | (_| | | \\__ \\  __/\n");
  printf("|___|_| |_|_|\\__|_|\\__,_|_|_|___/\\___|\n");
  printf("---------------------------------------\n");
}

//
// Print a divider to separate initialisation from computation
//
void print_integrate_divider(void) {
  printf("-------------------------------------------\n");
  printf(" ___       _                       _       \n");
  printf("|_ _|_ __ | |_ ___  __ _ _ __ __ _| |_ ___ \n");
  printf(" | || '_ \\| __/ _ \\/ _` | '__/ _` | __/ _ \\\n");
  printf(" | || | | | ||  __/ (_| | | | (_| | ||  __/\n");
  printf("|___|_| |_|\\__\\___|\\__, |_|  \\__,_|\\__\\___|\n");
  printf("                   |___/                    \n");
  printf("-------------------------------------------\n");
}

//
// Compute the average error in an array
//
double average_error(const double *expected, const double *actual, int size) {
  double sum_errors = 0.0;

  for (int i = 0; i < size; ++i) {
    double error = (expected[i] != 0.0) ? fabs((actual[i] - expected[i]) / expected[i]) : 0.0;
    sum_errors += error;
  }

  return sum_errors / size;
}

//
// Read in test temperatures from file
//
void load_temperatures(double **temperature_array, int *num_elements) {
  const char *filename = "test-temperatures.txt";
  FILE *file = fopen(filename, "r");
  if (!file) {
    perror("Error opening file");
    exit(EXIT_FAILURE);
  }

  // Count the number of lines in the file to determine the number of elements
  int count = 0;
  char c;
  while ((c = fgetc(file)) != EOF) {
    if (c == '\n') { count++; }
  }

  // Allocate memory for the array
  *temperature_array = malloc(count * sizeof(double));
  if ((*temperature_array) == NULL) {
    fclose(file);
    perror("Memory allocation failed");
    exit(EXIT_FAILURE);
  }

  // Reset file pointer to the beginning of the file
  fseek(file, 0, SEEK_SET);

  // Read temperatures into the array
  for (int i = 0; i < count; ++i) {
    if (fscanf(file, "%lf", &(*temperature_array)[i]) != 1) {
      fclose(file);
      free(*temperature_array);
      perror("Error reading temperature");
      exit(EXIT_FAILURE);
    }
  }

  fclose(file);
  *num_elements = count;
}

//
// Time how long it takes to compute alpha_sp for a given integrator function
//
double time_integrator(double (*integrator)(double (*integrand)(double, void *), void *, double, double, double),
                       double **results, int *results_count) {
  int num_temperatures;
  double *temperatures;

  load_temperatures(&temperatures, &num_temperatures);

  /* Count how big of an array we need to store the results for comparison */
  int count = 0;
  for (int i = 0; i < num_temperatures; ++i) {
    for (int j = 0; j < nlevels_macro; ++j) {
      for (int k = 0; k < xconfig[j].n_bfd_jump; ++k) { count++; }
    }
  }
  *results = calloc(count, sizeof(double));
  if ((*results) == NULL) {
    perror("Memory allocation failed");
    exit(EXIT_FAILURE);
  }
  if (results_count != NULL) { *results_count = count; }

  const clock_t start_time = clock();
  count = 0;
  for (int i = 0; i < num_temperatures; ++i) {
    const double temperature = temperatures[i];
    for (int j = 0; j < nlevels_macro; ++j) {
      for (int k = 0; k < xconfig[j].n_bfd_jump; ++k) {
        const double result = alpha_sp(&phot_top[xconfig[j].bfd_jump[k]], temperature, 0, integrator);
        (*results)[count] = result;
        count++;
      }
    }
  }
  const clock_t end_time = clock();
  free(temperatures);

  return ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
}

//
// Main function of the program
//
int main(int argc, char **argv) {

  geo.ioniz_mode = 9;
  print_initialise_divider();
  Log_set_verbosity(SHOW_LOG);
  get_atomic_data("data/h10_hetop_standard80.dat");

  Log_set_verbosity(SHOW_ERROR);
  print_integrate_divider();
  gsl_set_error_handler_off();

  int results_count;
  double *results_default;
  double *results_cquad;
  double *results_qag;
  double *results_qags;
  double *results_romberg;
  double *results_trap;

  const double time_default = time_integrator(integrate_default, &results_default, &results_count);
  printf("Default integrator    : %f seconds\n", time_default);

  const double time_trap = time_integrator(integrate_trap, &results_trap, NULL);
  printf("Trapezium integrator  : %f seconds : average fractional error %g\n", time_trap,
         average_error(results_default, results_trap, results_count));

  const double time_cquad = time_integrator(integrate_cquad, &results_cquad, NULL);
  printf("CQUAD integrator      : %f seconds : average fractional error %g\n", time_cquad,
         average_error(results_default, results_cquad, results_count));

  const double time_qag = time_integrator(integrate_qag, &results_qag, NULL);
  printf("QAG integrator        : %f seconds : average fractional error %g\n", time_qag,
         average_error(results_default, results_qag, results_count));

  const double time_qags = time_integrator(integrate_qags_small, &results_qags, NULL);
  printf("Small QAGS integrator : %f seconds : average fractional error %g\n", time_qags,
         average_error(results_default, results_qags, results_count));

  const double time_romberg = time_integrator(integrate_romberg, &results_romberg, NULL);
  printf("Romberg integrator    : %f seconds : average fractional error %g\n", time_romberg,
         average_error(results_default, results_romberg, results_count));

  free(results_default);
  free(results_cquad);
  free(results_qag);
  free(results_romberg);

  return EXIT_SUCCESS;
}
