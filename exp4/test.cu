#include <stdio.h>

__global__ void mykernel()
{
    printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
}

int main()
{
    mykernel<<<2, 3>>>();
    cudaDeviceSynchronize();
    return 0;
}