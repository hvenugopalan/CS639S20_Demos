all: files
files: main.cpp Laplacian.cpp Reductions.cpp PointwiseOps.cpp Utilities.cpp MatVecMultiply.cpp Substitutions.cpp ConjugateGradients.cpp *.h
	g++ -o LaplacianSolver main.cpp Laplacian.cpp Reductions.cpp PointwiseOps.cpp Utilities.cpp MatVecMultiply.cpp Substitutions.cpp ConjugateGradients.cpp *.h -m64 -I/u/h/a/harikrishnan/intel/mkl/include/ -L/u/h/a/harikrishnan/intel/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -fopenmp
clean:
	rm -f *.o LaplacianSolver
	rm -f *.pgm
