#include "Saxpy_Norm.h"
#include <algorithm>

float Saxpy_Norm(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const float scale)
{
    float result = 0.; 

    #pragma omp parallel for reduction(max:result)
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
    {
        z[i][j][k] = x[i][j][k] * scale + y[i][j][k];
        result = std::max(result, std::abs(z[i][j][k]));

    }

    return result;
}
            