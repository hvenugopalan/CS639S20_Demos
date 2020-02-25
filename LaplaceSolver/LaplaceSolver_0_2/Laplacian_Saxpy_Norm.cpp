#include "Laplacian_Saxpy_Norm.h"
#include <algorithm>

float ComputeLaplacian_Saxpy_Norm(const float (&u)[XDIM][YDIM][ZDIM], 
                                float (&Lu)[XDIM][YDIM][ZDIM],
                                const float (&y)[XDIM][YDIM][ZDIM],
                                float (&z)[XDIM][YDIM][ZDIM],
                                const float scale)
{  
    float result = 0.;

    #pragma omp parallel for reduction(max:result)
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
    {
        Lu[i][j][k] =
            -6 * u[i][j][k]
            + u[i+1][j][k]
            + u[i-1][j][k]
            + u[i][j+1][k]
            + u[i][j-1][k]
            + u[i][j][k+1]
            + u[i][j][k-1];
        
        z[i][j][k] = Lu[i][j][k] * scale + y[i][j][k];
        result = std::max(result, std::abs(z[i][j][k]));

    }
    
    return result;
            
}
