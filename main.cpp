#include "ConjugateGradients.h"
#include "Laplacian.h"
#include "Timer.h"
#include "Utilities.h"

Timer timerLaplacian;
Timer timerSaxpy1;
Timer timerSaxpy2;
Timer timerNorm;
Timer timerInnerProduct;
Timer timerCopy;

int main(int argc, char *argv[])
{
    using array_t = float (&) [XDIM][YDIM][ZDIM];

    float *xRaw = new float [XDIM*YDIM*ZDIM];
    float *fRaw = new float [XDIM*YDIM*ZDIM];
    float *pRaw = new float [XDIM*YDIM*ZDIM];
    float *rRaw = new float [XDIM*YDIM*ZDIM];
    float *zRaw = new float [XDIM*YDIM*ZDIM];
    
    array_t x = reinterpret_cast<array_t>(*xRaw);
    array_t f = reinterpret_cast<array_t>(*fRaw);
    array_t p = reinterpret_cast<array_t>(*pRaw);
    array_t r = reinterpret_cast<array_t>(*rRaw);
    array_t z = reinterpret_cast<array_t>(*zRaw);
    
    CSRMatrix matrix;
    CSRMatrix L;

    // Initialization
    {
        Timer timer;
        timer.Start();
        InitializeProblem(x, f);
        
        L = BuildPreconditionerMatrix(); // This takes a while ...

        // Make sure that the preconditioner is a valid lower-triangular CSR matrix ...
        L.CheckLowerTriangular();
        
        timer.Stop("Initialization : ");
    }

    // Call Conjugate Gradients algorithm
    {	
    timerLaplacian.Reset(); timerSaxpy1.Reset(); timerSaxpy2.Reset();	timerNorm.Reset(); timerInnerProduct.Reset(); timerCopy.Reset();
    ConjugateGradients( L, x, f, p, r, z, true);
    timerLaplacian.Print("Total Laplacian Time : ");
    timerSaxpy1.Print("Total Saxpy1 Time : ");
    timerSaxpy2.Print("Total Saxpy2 Time : ");
    timerNorm.Print("Total Norm Time : ");
    timerInnerProduct.Print("Total Inner Product Time : ");
    timerCopy.Print("Total Copy Time : ");
    }

    return 0;
}