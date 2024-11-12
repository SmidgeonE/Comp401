#include <iostream>
#include <cudaExample.cuh>

__global__ void addThreeNumbers(int a, int b, int c, int *result) {
    *result = a + b + c;
}

int main() {
    int a = 5, b = 10, c = 15;
    int result;
    int *d_result;

    cudaMalloc((void**)&d_result, sizeof(int));

    addThreeNumbers<<<1, 1>>>(a, b, c, d_result);

    cudaMemcpy(&result, d_result, sizeof(int), cudaMemcpyDeviceToHost);

    std::cout << "The result of adding " << a << ", " << b << ", and " << c << " is " << result << std::endl;

    cudaFree(d_result);

    return 0;
}