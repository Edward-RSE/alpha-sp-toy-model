#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <atomic.h>
#include <python.h>

int main(void) {
  MPI_Init(NULL, NULL);

  int host_rank;
  int host_size;
  MPI_Comm host_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &host_comm);
  MPI_Comm_rank(host_comm, &host_rank);
  MPI_Comm_size(host_comm, &host_size);

  //
  // Example with an array of 5 integers
  //

  int root_data[5] = {1, 2, 3, 4, 5};
  int disp_unit = (int) sizeof(int);
  MPI_Aint window_size = 0;
  if (host_rank == 0) {
    window_size = 5 * disp_unit;// sharing a single piece of data
  }

  MPI_Win win;
  int *shared_data;
  MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, host_comm, &shared_data, &win);

  if (host_rank != 0) { MPI_Win_shared_query(win, 0, &window_size, &disp_unit, &shared_data); }

  if (host_rank == 0) {
    for (int i = 0; i < 5; ++i) { shared_data[i] = root_data[i]; }
  }

  MPI_Win_fence(0, win);

  for (int i = 0; i < host_size; ++i) {
    if (i == host_rank) {
      printf("Host rank %d has values [", host_rank);
      for (int j = 0; j < 5; ++j) { printf(" %d", shared_data[j]); }
      printf(" ]\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Win_free(&win);

  //
  // Example with Python elements struct
  //

  Log_set_verbosity(0);
  if (host_rank == 0) { get_atomic_data("data/standard80.dat"); }

  struct elements *shared_ele;
  disp_unit = (int) sizeof(ele_dummy);
  window_size = 0;

  if (host_rank == 0) {
    window_size = nelements * disp_unit;
    MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, host_comm, &shared_ele, &win);
    memcpy(shared_ele, ele, window_size);
    free(ele);
  } else {
    MPI_Win_allocate_shared(0, disp_unit, MPI_INFO_NULL, host_comm, &shared_ele, &win);
    MPI_Win_shared_query(win, 0, &window_size, &disp_unit, &shared_ele);
  }

  MPI_Win_fence(0, win);
  for (int i = 0; i < host_size; ++i) {
    if (i == host_rank) {
      printf("Host rank %d: element %d name = %s\n", host_rank, host_rank, shared_ele[host_rank].name);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Win_fence(0, win);

  MPI_Win_free(&win);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
