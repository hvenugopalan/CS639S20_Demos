#include "Combined_Saxpy.h"

void Saxpy(float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM],
    const float (&z)[XDIM][YDIM][ZDIM], const float alpha, const float beta)
{
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
    {
        y[i][j][k] = x[i][j][k] * alpha + y[i][j][k];
        x[i][j][k] = x[i][j][k] * beta + z[i][j][k];
    }
}
