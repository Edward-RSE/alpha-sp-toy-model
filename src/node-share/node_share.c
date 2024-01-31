#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <atomic.h>
#include <python.h>

int main(void) {
  MPI_Init(NULL, NULL);

  // In this toy model, we'll initialise memory in a root rank and share it with
  // all other non-root ranks. This only works with intra-node processes, e.g.
  // processes on the same machine. So data can be shared with this method with
  // all CPUs on the same machine (even those in different NUMA regions), but
  // not with other machines on the same network (e.g. other nodes in a HPC
  // cluster). Therefore, we can no longer really use MPI_COMM_WORLD and have to
  // split up MPI into other communicators where we have one communicator PER
  // NODE. The MPI processes running on ONE NODE will belong to the same
  // communicator. We can split MPI_COMM_WORLD into something like this
  // automatically using MPI_Comm_split_type with the MPI_COMM_TYPE_SHARED
  // option.
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

  // Variables for shared window buffer stuff
  MPI_Win win;
  int disp_unit = (int) sizeof(int);
  MPI_Aint window_size = 0;
  int *shared_data;  // The pointer to shared data

  // On the root rank for the host, we'll create a shared memory buffer of the
  // appropriate size in `shared_data`. This will also allocate space for the
  // `win` handle. On non-root ranks, we do not allocate space for shared data
  // but we still need to set up the window. Using `MPI_Win_shared_query` we can
  // point `shared_data` to the same address as the shared buffer on the root
  // rank to gain access to the data.
  if (host_rank == 0) {
    window_size = 5 * disp_unit;
    MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, host_comm, &shared_data, &win);
  } else {
    MPI_Win_allocate_shared(0, disp_unit, MPI_INFO_NULL, host_comm, &shared_data, &win);
    MPI_Win_shared_query(win, 0, &window_size, &disp_unit, &shared_data);
  }
  MPI_Win_fence(0, win);  // Force synchronisation to make sure all windows are in the same state

  // At this point, shared_data will be (hopefully) full of zero's. But we can
  // modify the contents of it and this should be reflected in the pointer on
  // non-root ranks. This is not to say it's making a copy in each rank, but
  // because the non-root ranks have a pointer to memory in the root rank then
  // they can see the modifications made to the (shared) data in the root rank
  if (host_rank == 0) {
    for (int i = 0; i < 5; ++i) { shared_data[i] = root_data[i]; }
  }
  MPI_Win_fence(0, win);  // Again, force synchronisation

  // Now when we come to print the contents of shared_data, each rank should
  // have the same data even though we have only explicitly set the value of
  // the data in the root rank
  for (int i = 0; i < host_size; ++i) {
    if (i == host_rank) {
      printf("Host rank %d has values [", host_rank);
      for (int j = 0; j < 5; ++j) { printf(" %d", shared_data[j]); }
      printf(" ]\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Force a synchronisation and free the window after we have done with the
  // data in this shared window
  MPI_Win_fence(0, win);
  MPI_Win_free(&win);

  //
  // Example with Python elements struct
  //

  Log_set_verbosity(SHOW_LOG);

  // Only read in data for the host rank
  if (host_rank == 0) { get_atomic_data("data/standard80.dat"); }

  // Create variables for window share
  struct elements *shared_ele;
  disp_unit = (int) sizeof(ele_dummy);
  window_size = 0;

  // On the root host rank, allocate shared memory space for shared_ele of the
  // appropriate window size. Then copy the contents of ele into shared_ele and
  // free ele.
  // For non-root ranks, we create the shared pointer and then used shared_query
  // to point it to the correct base address in the root rank
  if (host_rank == 0) {
    window_size = nelements * disp_unit;
    MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, host_comm, &shared_ele, &win);
    memcpy(shared_ele, ele, window_size);
    free(ele);
  } else {
    MPI_Win_allocate_shared(0, disp_unit, MPI_INFO_NULL, host_comm, &shared_ele, &win);
    MPI_Win_shared_query(win, 0, &window_size, &disp_unit, &shared_ele);
  }
  MPI_Win_fence(0, win); // force synchronisation, to make sure all windows are the same

  // Point the shared_ele data back to ele
  ele = shared_ele;

  // Print stuff
  for (int i = 0; i < host_size; ++i) {
    if (i == host_rank) {
      printf("Host rank %d: element %d name = %s\n", host_rank, host_rank, ele[host_rank].name);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Synchronise and free the window at the end
  MPI_Win_fence(0, win);
  MPI_Win_free(&win);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
