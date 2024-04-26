#include <chrono>
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <cstring>
#include <cmath>
#include <algorithm>

#define EPS 1e-5

namespace ch = std::chrono;

void Ring_Allreduce(void *sendbuf, void *recvbuf, int n, MPI_Comm comm, int comm_sz, int my_rank)
{
    // TODO
    const int P = comm_sz;
    const int FULL_BLOCK_SIZE = (n + P - 1) / P;
    const int EDGE_BLOCK_SIZE = ((n + FULL_BLOCK_SIZE - 1) % FULL_BLOCK_SIZE) + 1;
    float *my_send_array = (float *)sendbuf;
    float *my_recv_array = (float *)recvbuf;

    memcpy(my_recv_array + my_rank * FULL_BLOCK_SIZE,
           my_send_array + my_rank * FULL_BLOCK_SIZE,
           my_rank == P - 1 ? (EDGE_BLOCK_SIZE * sizeof(float)) : (FULL_BLOCK_SIZE * sizeof(float)));

    for (int k = 0; k < P; ++k)
    {
        const int SEND_RANK = (my_rank - k + P) % P;
        const int SEND_SIZE = SEND_RANK == P - 1 ? EDGE_BLOCK_SIZE : FULL_BLOCK_SIZE;
        const int RECV_RANK = (((my_rank - 1 + P) % P) - k + P) % P;
        const int RECV_SIZE = RECV_RANK == P - 1 ? EDGE_BLOCK_SIZE : FULL_BLOCK_SIZE;
        MPI_Sendrecv(
            &my_recv_array[SEND_RANK * FULL_BLOCK_SIZE], SEND_SIZE, MPI_FLOAT, (my_rank + 1) % P, 0,
            &my_recv_array[RECV_RANK * FULL_BLOCK_SIZE], RECV_SIZE, MPI_FLOAT, (my_rank - 1 + P) % P, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < RECV_SIZE; ++i)
            my_recv_array[RECV_RANK * FULL_BLOCK_SIZE + i] += my_send_array[RECV_RANK * FULL_BLOCK_SIZE + i];
    }

    for (int k = 0; k < P; ++k)
    {
        const int SEND_RANK = (my_rank + 1 - k + P) % P;
        const int SEND_SIZE = SEND_RANK == P - 1 ? EDGE_BLOCK_SIZE : FULL_BLOCK_SIZE;
        const int RECV_RANK = (my_rank - k + P) % P;
        const int RECV_SIZE = RECV_RANK == P - 1 ? EDGE_BLOCK_SIZE : FULL_BLOCK_SIZE;
        MPI_Sendrecv(
            &my_recv_array[SEND_RANK * FULL_BLOCK_SIZE], SEND_SIZE, MPI_FLOAT, (my_rank + 1) % P, 0,
            &my_recv_array[RECV_RANK * FULL_BLOCK_SIZE], RECV_SIZE, MPI_FLOAT, (my_rank - 1 + P) % P, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

// reduce + bcast
void Naive_Allreduce(void *sendbuf, void *recvbuf, int n, MPI_Comm comm, int comm_sz, int my_rank)
{
    MPI_Reduce(sendbuf, recvbuf, n, MPI_FLOAT, MPI_SUM, 0, comm);
    MPI_Bcast(recvbuf, n, MPI_FLOAT, 0, comm);
}

int main(int argc, char *argv[])
{
    int ITER = atoi(argv[1]);
    int n = atoi(argv[2]);
    float *mpi_sendbuf = new float[n];
    float *mpi_recvbuf = new float[n];
    float *naive_sendbuf = new float[n];
    float *naive_recvbuf = new float[n];
    float *ring_sendbuf = new float[n];
    float *ring_recvbuf = new float[n];

    MPI_Init(nullptr, nullptr);
    int comm_sz;
    int my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    srand(time(NULL) + my_rank);
    for (int i = 0; i < n; ++i)
        mpi_sendbuf[i] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    memcpy(naive_sendbuf, mpi_sendbuf, n * sizeof(float));
    memcpy(ring_sendbuf, mpi_sendbuf, n * sizeof(float));

    // warmup and check
    MPI_Allreduce(mpi_sendbuf, mpi_recvbuf, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    Naive_Allreduce(naive_sendbuf, naive_recvbuf, n, MPI_COMM_WORLD, comm_sz, my_rank);
    Ring_Allreduce(ring_sendbuf, ring_recvbuf, n, MPI_COMM_WORLD, comm_sz, my_rank);
    bool correct = true;
    for (int i = 0; i < n; ++i)
        if (abs(mpi_recvbuf[i] - ring_recvbuf[i]) > EPS)
        {
            correct = false;
            break;
        }

    if (correct)
    {
        auto beg = ch::high_resolution_clock::now();
        for (int iter = 0; iter < ITER; ++iter)
            MPI_Allreduce(mpi_sendbuf, mpi_recvbuf, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        auto end = ch::high_resolution_clock::now();
        double mpi_dur = ch::duration_cast<ch::duration<double>>(end - beg).count() * 1000; // ms

        beg = ch::high_resolution_clock::now();
        for (int iter = 0; iter < ITER; ++iter)
            Naive_Allreduce(naive_sendbuf, naive_recvbuf, n, MPI_COMM_WORLD, comm_sz, my_rank);
        end = ch::high_resolution_clock::now();
        double naive_dur = ch::duration_cast<ch::duration<double>>(end - beg).count() * 1000; // ms

        beg = ch::high_resolution_clock::now();
        for (int iter = 0; iter < ITER; ++iter)
            Ring_Allreduce(ring_sendbuf, ring_recvbuf, n, MPI_COMM_WORLD, comm_sz, my_rank);
        end = ch::high_resolution_clock::now();
        double ring_dur = ch::duration_cast<ch::duration<double>>(end - beg).count() * 1000; // ms

        if (my_rank == 0)
        {
            std::cout << "Correct." << std::endl;
            std::cout << "MPI_Allreduce:   " << mpi_dur << " ms." << std::endl;
            std::cout << "Naive_Allreduce: " << naive_dur << " ms." << std::endl;
            std::cout << "Ring_Allreduce:  " << ring_dur << " ms." << std::endl;
        }
    }
    else if (my_rank == 0)
        std::cout << "Wrong!" << std::endl;

    delete[] mpi_sendbuf;
    delete[] mpi_recvbuf;
    delete[] naive_sendbuf;
    delete[] naive_recvbuf;
    delete[] ring_sendbuf;
    delete[] ring_recvbuf;
    MPI_Finalize();
    return 0;
}