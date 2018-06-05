//     #include <functions.h>

// Headers
#include <complex.h> // Manipulating compelx numbers
#include <float.h>   // Various limits (float and double)
#include <limits.h>  // Various limits (integer)
#include <math.h>    // Common mathematical functions
#include <stdio.h>   // Core input/output operations
#include <stdlib.h>  // Conversions, random numbers, memory allocation, etc.
#include <string.h>  // String handling functions
#include <time.h>    // Converting between various date/time formats
#include <unistd.h>  // Hostname, etc.

// Conditional inclusions
// OpenMP parallelization
#ifdef PAR_OMP
  #include <omp.h>
#endif

// MPICH parallelization
#ifdef PAR_MPI 
  #include <mpi.h>
#endif

// Hybrid MPICH+OpenMP parallelization
#ifdef PAR_HYB
  #include <omp.h>
  #include <mpi.h>
#endif

