#include <mpi.h>
#include <stdio.h>

double f(double x)
{
    double y = 1 + x + x * x;
    return y;
}

double Trap(double left_endpt, double right_endpt, int trap_count, double base_len)
{
    double estimate, x;
    int i;

    estimate = (f(left_endpt) + f(right_endpt)) / 2.0;
    for (i = 1; i <= trap_count - 1; i++)
    {
        x = left_endpt + i * base_len;
        estimate += f(x);
    }
    estimate = estimate * base_len;

    return estimate;
} /** Trap */

int main(int argc, char **argv)
{
    int my_rank, comm_size, n = 1024, local_n;
    double a = 0.0, b = 3.0, h, local_a, local_b;
    double local_integral, total_integral;
    int source;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    h = (b - a) / n;
    local_n = n / comm_size;
    local_a = a + my_rank * local_n * h;
    local_b = local_a + local_n * h;

    local_integral = Trap(local_a, local_b, local_n, h);

    // if (my_rank != 0)
    // {
    //     MPI_Send(&local_integral, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    // }
    // else
    // {
    //     total_integral = local_integral;
    //     for (source = 1; source < comm_size; source++)
    //     {
    //         MPI_Recv(&local_integral, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         total_integral += local_integral;
    //     }
    // }

    total_integral = 0;
    MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        printf("With n = %d trapezoids, our estimate\n", n);
        printf("of the intergral from %f to %f = %.15e\n", a, b, total_integral);
    }

    MPI_Finalize();
    return 0;
}