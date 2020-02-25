#pragma once

#include "Parameters.h"


float ComputeLaplacian_InnerProduct(const float (&u)[XDIM][YDIM][ZDIM], 
                                float (&Lu)[XDIM][YDIM][ZDIM]);
