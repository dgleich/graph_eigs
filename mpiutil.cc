/**
 * @file mpiutil.cc
 * ================
 * Handy utilities to work with MPI.
 *
 * @author David F. Gleich
 */
 
/**
 * History
 * ------- 
 * :2010-11-07: Based on mcrapr functions from MPI monte-carlo experiments.
 */
 
#include <stdio.h>
#include <mpi.h>
 
int mpi_printf(MPI_Comm comm, const char format[], ...) 
{
    int rank;
    MPI_Comm_rank(comm,&rank); 
    if (!rank) {
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
        fflush(stdout);
    }
    return (0);
}

int mpi_world_printf(const char format[], ...)
{   
    // TODO figure out how to have this one call mpi_printf
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    if (!rank) {
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
        fflush(stdout);
    }
    return (0);
}

int mpi_sync_vprintf(MPI_Comm comm, const char format[], va_list args)
{
    int rank, size;
    int token=0;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);
    
    MPI_Barrier(comm);
    
    if (rank) {
        MPI_Status s;
        MPI_Recv(&token,1,MPI_INT,rank-1,0,comm,&s);
    }

    // always execute this code
    vprintf(format, args);
    fflush(stdout);
    
    if (rank+1!=size) {
        MPI_Send(&token,1,MPI_INT,rank+1,0,comm);
    }

    MPI_Barrier(comm);
    
    return (0);
}


 
/** A synchronized printf prints from each node in sequence.
 */
int mpi_sync_printf(MPI_Comm comm, const char format[], ...) 
{
    // always execute this code
    va_list args;
    va_start(args, format);
    mpi_sync_vprintf(comm, format, args);
    va_end(args);
    
    return (0);
} 

/** A synchronized printf prints from each node in sequence.
 */
int mpi_world_sync_printf(const char format[], ...) 
{
    // always execute this code
    va_list args;
    va_start(args, format);
    mpi_sync_vprintf(MPI_COMM_WORLD, format, args);
    va_end(args);
    
    return (0);
}
