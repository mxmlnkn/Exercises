#include "mpi.h"
#include <stdio.h>
#include <string.h>
 
int main(int argc, char *argv[])
{
    int myrank, message_size=50, tag=99;
    char message[message_size];
    MPI_Status status;
 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
    if (myrank == 0) {
        MPI_Recv(message, message_size, MPI_CHAR, 1, tag, MPI_COMM_WORLD, &status);
        printf("received \"%s\"\n", message);
    }
    else {
        strcpy(message, "Hello, there");
        MPI_Send(message, strlen(message)+1, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}