#pragma once

#include "Parameters.h"

// Scale array x by given number, add y, write result into z and take norm of z
float Saxpy_Norm(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM], const float scale);
