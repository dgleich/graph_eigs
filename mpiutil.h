#ifndef MPIUTIL_H
#define MPIUTIL_H

#include <stdio.h>
 
int mpi_printf(MPI_Comm comm, const char format[], ...);
int mpi_world_printf(const char format[], ...);

int mpi_sync_vprintf(MPI_Comm comm, const char format[], va_list args);
int mpi_sync_printf(MPI_Comm comm, const char format[], ...);
int mpi_world_sync_printf(const char format[], ...);

#endif /* MPIUTIL_H */
