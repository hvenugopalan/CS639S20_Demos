#pragma once

#include "Parameters.h"

// Scale array x by alpha, add y, and write result into y
// Scale array x by beta, add z, and write result into x
void Saxpy(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM], const float alpha, const float beta);
