#pragma once

#include "Parameters.h"


float ComputeLaplacian_Saxpy_Norm(const float (&u)[XDIM][YDIM][ZDIM], 
                                float (&Lu)[XDIM][YDIM][ZDIM],
                                const float (&y)[XDIM][YDIM][ZDIM],
                                float (&z)[XDIM][YDIM][ZDIM],
                                const float scale);
