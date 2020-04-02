#include "ConjugateGradients.h"
#include "Laplacian.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Substitutions.h"
#include "Utilities.h"
#include "Timer.h"

#include <iostream>

extern Timer timerLaplacian;
extern Timer timerSaxpy1;
extern Timer timerSaxpy2;
extern Timer timerNorm;
extern Timer timerInnerProduct;
extern Timer timerCopy;

void ConjugateGradients(
    CSRMatrix& L,
    float (&x)[XDIM][YDIM][ZDIM],
    const float (&f)[XDIM][YDIM][ZDIM],
    float (&p)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const bool writeIterations)
{
    // Algorithm : Line 2
    timerLaplacian.Restart(); ComputeLaplacian( x, z); timerLaplacian.Pause();
    timerSaxpy2.Restart(); Saxpy(z, f, r, -1); timerSaxpy2.Pause();
    float nu = Norm(r);

    // Algorithm : Line 3
    if (nu < nuMax) return;

    // Algorithm : Line 4
    timerCopy.Restart(); Copy(r, p); timerCopy.Pause();
    ForwardSubstitution(L, &p[0][0][0]);
    BackwardSubstitution(L, &p[0][0][0]);
    timerInnerProduct.Restart();
    float rho=InnerProduct(p, r);
    timerInnerProduct.Pause();
        
    // Beginning of loop from Line 5
    for(int k=0;;k++)
    {
        std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
        timerLaplacian.Restart(); ComputeLaplacian(p, z); timerLaplacian.Pause();
        timerInnerProduct.Restart();
	float sigma=InnerProduct(p, z);
	timerInnerProduct.Pause();

        // Algorithm : Line 7
        float alpha=rho/sigma;

        // Algorithm : Line 8
        timerSaxpy1.Restart(); Saxpy(z, r, -alpha); timerSaxpy1.Pause();
        timerNorm.Restart();
	nu=Norm(r);
	timerNorm.Pause();

        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
            timerSaxpy1.Restart(); Saxpy(p, x, alpha); timerSaxpy1.Pause();
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            return;
        }

        // Algorithm : Line 13
	timerCopy.Restart();
        Copy(r, z);
	timerCopy.Pause();
        ForwardSubstitution(L, &z[0][0][0]);
        BackwardSubstitution(L, &z[0][0][0]);
        timerInnerProduct.Restart();
	float rho_new = InnerProduct(z, r);
	timerInnerProduct.Pause();

        // Algorithm : Line 14
        float beta = rho_new/rho;

        // Algorithm : Line 15
        rho=rho_new;

        // Algorithm : Line 16
        timerSaxpy1.Restart(); Saxpy(p, x, alpha); timerSaxpy1.Pause();
        timerSaxpy2.Restart(); Saxpy(p, z, p, beta); timerSaxpy2.Pause();

        if (writeIterations) WriteAsImage("x", x, k, 0, XDIM/2);
    }

}
