#include <mpi.h>
#include <stdio.h>
#include <string.h>

#define MAX_STRING 100

int main(int argc, char **argv)
{
    char greeting[MAX_STRING];
    int comm_sz, my_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    if (my_rank != 0)
    {
        sprintf(greeting, "Greentings from process %d of %d!", my_rank, comm_sz);
        MPI_Send(greeting, strlen(greeting) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        printf("Greentings from process %d of %d!\n", my_rank, comm_sz);

        for (int s = 1; s < comm_sz; s++)
        {
            MPI_Recv(greeting, MAX_STRING, MPI_CHAR, s, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("%s\n",greeting);
        }
    }

    MPI_Finalize();
    return 0;
}