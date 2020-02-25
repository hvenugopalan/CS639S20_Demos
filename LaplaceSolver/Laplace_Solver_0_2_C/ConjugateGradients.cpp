#include "Laplacian.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Utilities.h"
#include "Timer.h"
#include <iostream>

void ConjugateGradients(
    float (&x)[XDIM][YDIM][ZDIM],
    const float (&f)[XDIM][YDIM][ZDIM],
    float (&p)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const bool writeIterations)
{
    Timer timer_cl_2,timer_cl_6;
    Timer timer_s_2, timer_s_8, timer_s_9, timer_s_16_1, timer_s_16_2;
    Timer timer_n_2, timer_n_8;
    Timer timer_ip_4, timer_ip_6, timer_ip_13;
    Timer timer_cg;
    // Algorithm : Line 2

    timer_cg.Start();
    timer_cl_2.Start();
    ComputeLaplacian(x, z);
    timer_cl_2.Pause();
    
    timer_s_2.Restart();
    Saxpy(z, f, r, -1);
    timer_s_2.Pause();
    
    timer_n_2.Restart();
    float nu = Norm(r);
    timer_n_2.Pause();

    // Algorithm : Line 3
    if (nu < nuMax) return;
        
    // Algorithm : Line 4
    Copy(r, p);
    timer_ip_4.Restart();
    float rho=InnerProduct(p, r);
    timer_ip_4.Pause();
        
    // Beginning of loop from Line 5
    for(int k=0;;k++)
    {
        //std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
        timer_cl_6.Restart();
        ComputeLaplacian(p, z);
        timer_cl_6.Pause();

        timer_ip_6.Restart();
        float sigma=InnerProduct(p, z);
        timer_ip_6.Pause();


        // Algorithm : Line 7
        float alpha=rho/sigma;

        // Algorithm : Line 8
        timer_s_8.Restart();
        Saxpy(z, r, r, -alpha);
        timer_s_8.Pause();

        timer_n_8.Restart();
        nu=Norm(r);
        timer_n_8.Pause();

        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
            timer_s_9.Restart();
            Saxpy(p, x, x, alpha);
            timer_s_9.Pause();

            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            timer_cl_2.Print("Compute Laplacian on line 2:");
            timer_cl_6.Print("Compute Laplacian on line 6:");
            timer_s_2.Print("Saxpy on line 2:");
            timer_s_8.Print("Saxpy on line 8:");
            timer_s_9.Print("Saxpy on line 9:");
            timer_s_16_1.Print("1st Saxpy on line 16:");
            timer_s_16_2.Print("2nd Saxpy on line 16:");
            timer_n_2.Print("Norm on line 2:");
            timer_n_8.Print("Norm on line 8:");
            timer_ip_4.Print("Inner Product on line 4:");
            timer_ip_6.Print("Inner Product on line 6:");
            timer_ip_13.Print("Inner Product on line 13:");
            timer_cg.Stop("ConjugateGradients: ");
            return;
        }
            
        // Algorithm : Line 13
        Copy(r, z);
        timer_ip_13.Restart();
        float rho_new = InnerProduct(z, r);
        timer_ip_13.Pause();

        // Algorithm : Line 14
        float beta = rho_new/rho;

        // Algorithm : Line 15
        rho=rho_new;

        // Algorithm : Line 16
        timer_s_16_1.Restart();
        Saxpy(p, x, x, alpha);
        timer_s_16_1.Pause();

        timer_s_16_2.Restart();
        Saxpy(p, r, p, beta);
        timer_s_16_2.Pause();

        if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
}
