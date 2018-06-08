#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include <mkl_lapacke.h>

#ifdef MODE_DEBUG
  #define N 20
  #define Ni 10
  #define Ne 10
  #define num_node 5

 /* Parameters */
  #define NN 3
  #define NRHS 1
  #define LDA NN
  #define LDB NN
#endif

#ifdef MODE_PROD
  #define N 1000
  #define Ni 500
  #define Ne 500
  #define num_node 100
  
 /* Parameters */
  #define NN 98
  #define NRHS 1
  #define LDA NN
  #define LDB NN
#endif

#define PI 3.141592653589793

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT mm, MKL_INT nn, double* aa, MKL_INT lda );
extern void print_int_vector( char* desc, MKL_INT nn, MKL_INT* aa );

// define the user-defined functions

// Electric Potential
double get_phi(double qn[]);

// Electric Field
double get_field_pic(double phi[],double dn);

// Linear Interpolation 
double interp_linear( double X[], double E[], double dn,double n[] );

double E[N];
double Ex[N]; 
double bb[NN];

int main(int argc, char *argv[])
{
  int i, j;
  double t, a1 = 0.0, b1 = 1.0;
  
  // Memory Preallocation
  double *node   = (double *) malloc((N) * sizeof(double));
  if (node == NULL) {
  printf("  Not enough memory to assign node[%d]\n", N);
  exit(1);
  }else {
  printf("\n");
  printf("  Enough memory to assign node[%d]\n", N);
  double node[N];
  // memset( node, 0.0, N * sizeof(double) );
  }
  double *xi   = (double *) malloc((N) * sizeof(double));
  if (xi == NULL) {
  printf("  Not enough memory to assign xi[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign xi[%d]\n", N);
  double xi[N];
  // memset( xi, 0.0, N * sizeof(double) );
  }
  double *xe   = (double *) malloc((N) * sizeof(double));
  if (xe == NULL) {
  printf("  Not enough memory to assign xe[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign xe[%d]\n", N);
  double xe[N];
  // memset( xe, 0.0, N * sizeof(double) );
  }
  double *vi   = (double *) malloc((N) * sizeof(double));
  if (vi == NULL) {
  printf("  Not enough memory to assign vi[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign vi[%d]\n", N);
  double vi[N];
  // memset( vi, 0.0, N * sizeof(double) );
  }
  double *ve   = (double *) malloc((N) * sizeof(double));
  if (ve == NULL) {
  printf("  Not enough memory to assign ve[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign ve[%d]\n", N);
  double ve[N];
  // memset( ve, 0.0, N * sizeof(double) );
  }
  double *X   = (double *) malloc((N) * sizeof(double));
  if (X == NULL) {
  printf("  Not enough memory to assign X[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign X[%d]\n", N);
  double X[N];
  // memset( X, 0.0, N * sizeof(double) );
  }
  double *V   = (double *) malloc((N) * sizeof(double));
  if (V == NULL) {
  printf("  Not enough memory to assign V[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign V[%d]\n", N);
  double V[N];
  // memset( V, 0.0, N * sizeof(double) );
  }
  double *Vnew   = (double *) malloc((N) * sizeof(double));
  if (Vnew == NULL) {
  printf("  Not enough memory to assign Vnew[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign Vnew[%d]\n", N);
  double Vnew[N];
  // memset( Vnew, 0.0, N * sizeof(double) );
  }
  double *Xnew   = (double *) malloc((N) * sizeof(double));
  if (Xnew == NULL) {
  printf("  Not enough memory to assign Xnew[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign Xnew[%d]\n", N);
  double Xnew[N];
  // memset( Xnew, 0.0, N * sizeof(double) );
  }
  //double *hx   = (double *) malloc((N) * sizeof(double));
  //if (hx == NULL) {
  //printf("  Not enough memory to assign hx[%d]\n", N);
  //exit(1);
  //}else { 
  //printf("\n");
  //printf("  Enough memory to assign hx[%d]\n", N);
  //double hx[N];
  //// memset( hx, 0.0, N * sizeof(double) );
  //}
  double hx[N][2];
  double w[N][2];
  for(i = 0;i<N;i++)
    {
    for(j=0;j<2;j++)
      {
      hx[i][j]  = 0.0;
      w[i][j]   = 0.0;
      }

    }
  double *n   = (double *) malloc((N) * sizeof(double));
  if (n == NULL) {
  printf("  Not enough memory to assign n[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign n[%d]\n", N);
  double n[N];
  // memset( n, 0.0, N * sizeof(double) );
  }
 // double *w   = (double *) malloc((N) * sizeof(double));
 // if (w == NULL) {
 // printf("  Not enough memory to assign w[%d]\n", N);
 // exit(1);
 // }else { 
 // printf("\n");
 // printf("  Enough memory to assign w[%d]\n", N);
 // double w[N];
 // // memset( w, 0.0, N * sizeof(double) );
 // }
  double *q   = (double *) malloc((N) * sizeof(double));
  if (q == NULL) {
  printf("  Not enough memory to assign q[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign q[%d]\n", N);
  double q[N];
  // memset( q, 0.0, N * sizeof(double) );
  }
  double *qn   = (double *) malloc((N) * sizeof(double));
  if (qn == NULL) {
  printf("  Not enough memory to assign qn[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign qn[%d]\n", N);
  double qn[N];
  // memset( qn, 0.0, N * sizeof(double) );
  }
  double *phi   = (double *) malloc((N) * sizeof(double));
  if (phi == NULL) {
  printf("  Not enough memory to assign phi[%d]\n", N);
  exit(1);
  }else { 
  printf("\n");
  printf("  Enough memory to assign phi[%d]\n", N);
  double phi[N];
  // memset( phi, 0.0, N * sizeof(double) );
  }

  for(i = 0; i < N; i++)
    {
    xi[i]   = 0.0;
    xe[i]   = 0.0;
    vi[i]   = 0.0;
    ve[i]   = 0.0;
    X[i]    = 0.0;
    V[i]    = 0.0;
    Xnew[i] = 0.0;    
    Vnew[i] = 0.0;
    //hx[i]   = 0.0;
    n[i]    = 0.0;
    //w[i]    = 0.0;
    q[i]    = 0.0;
    E[i]    = 0.0;
    Ex[i]   = 0.0;
    }
  for(i = 0; i < num_node; i++)
    {
    node[i] = 0.0;
    phi[i]  = 0.0;
    qn[i]    = 0.0;
    }

  double dxe, dxi, m = 1.0, k = 0.0, mi = 1000.0, me = 1.0;
  
  double dt = 0.0005;
  // double dt = 0.5;
  double r = 2;
  
  dxe = 2.0/Ne;
  dxi = 1.0/Ni; 

  // Initial Postion and Velocity of ions 
    for(i = 0 ; i< Ni; i++)
        {
            xi[i] = (m-0.5)*dxi;
            vi[i] = 0;
            m++;
        }

    m = 1.0;

    // Initial Postion and Velocity of electrons
    for( i = Ne; i < 3*Ne/2; i++)
        {
        xe[i] = (m-0.5)*dxe;
        ve[i] = 0.5 + 0.1 * sin(2.0*PI*xe[i]);
        m++;
        }
    for( i = 3*Ne/2; i < 2*Ne; i++)
        {
        xe[i] = xe[i-Ne/2];
        ve[i] = -0.5 - 0.1*sin(2.0*PI*xe[i]);
        }
    for(i = 0; i < N; i++)
        {
        X[i] = xi[i] + xe[i];
        V[i] = vi[i] + ve[i];
        #ifdef MODE_DEBUG
        printf("X = % 017.15lf   V = % 017.15lf\n",X[i],V[i]);
        #endif
        }

    // Creating 1D Mesh   
    for(i = 0; i < num_node; i++)
        {
        node[i] = k;
        #ifdef MODE_DEBUG
        printf("node(%d) = %.15lf\n",i,node[i]);
        #endif
        k = k + 1.0/(num_node-1.0);
        }
    
    double dn = 1.0/(num_node-1.0);

    // Time Loop BEGINS
    for(t = 0.0; t <= 1.5; t = t + dt) {

      printf("t = %6.5lf\n",t);

        // Calculating hx and Weight fraction
        for(i = 0;i<N;i++)
          {
          //n[i] = ceil(X[i]/dn);
          n[i] = floor((X[i] - a1)/(b1-a1)*(num_node - 1)+1);
          hx[i][1] =  X[i] - (n[i]-1.0)*dn;
          hx[i][2] = n[i]*dn - X[i];
          #ifdef MODE_DEBUG
          printf("n = %.15lf hx = %.15lf  %.15lf\n",n[i],hx[i][1],hx[i][2]);
          #endif
          }
        for(i = 1;i<N;i++)
          {
          w[i][1] = hx[i][1]/(hx[i][1]+hx[i][2]);
          w[i][2] = 1.0-w[i][1];
          #ifdef MODE_DEBUG
          printf("w1 = %.15lf  w2 = %.15lf\n",w[i][1],w[i][2]);
          #endif
          }
        // Random Charge generator
        srand((unsigned int)time(NULL));
        for(i=0;i<N/2;i++)
          {      
          q[i] = ((double)rand()/(double)(RAND_MAX)) * r;
          }
        for(i=N/2;i<N;i++)
          {
          q[i] = ((double)rand()/(double)(RAND_MAX)) * (-r);
          }
        #ifdef MODE_DEBUG
        for(i = 1;i<N;i++)
          {
          printf("q = %.15lf \n",q[i]);
          }
        #endif  
        // Calculating Charge density
        for(j = 1;j<N;j++)
          {
          n[j] = ceil(X[j]/dn);
          qn[(int)n[j]] = qn[(int)n[j]] +  w[j][2]*q[j];
          qn[(int)n[j]+1] = qn[(int)n[j]+1] +  w[j][1]*q[j];
          }
        #ifdef MODE_DEBUG  
          for(i=0;i<num_node;i++)
            {
            printf("qn = %.15lf \n",qn[i]);
            }
        #endif     
        // Calling Electric Potential Function
        get_phi(qn);
        #ifdef MODE_DEBUG
          for( i = 0; i < NN; i++ ) 
            {
            for( j = 0; j < NRHS; j++ ) 
              {
              printf( " %6.2lf", bb[i+j*LDA] );
              }
            }
        #endif

        for(i=1;i<num_node-1;i++)
          {
          phi[i] = bb[i-1];
          }
        // Calling Electric Field Function
        get_field_pic(phi,dn);
        #ifdef MODE_DEBUG
          for(i = 0;i < num_node; i++)
            {
            printf("phi = %.15lf\t E = %.15lf \n",phi[i],E[i]);
            }
        #endif       
        // Calling Linear Interpolation Function
        interp_linear( X, E, dn, n );
        #ifdef MODE_DEBUG
          for(i = 0; i < N; i++)
            {
            printf("Ex = %.15lf \n",Ex[i]);
            }
        #endif       
        
        // Updating Velocity
        for(i = 0;i < N; i++)
          {
          if(q[i] > 0)
              Vnew[i] = V[i] + dt * q[i] * Ex[i]/mi;
          else
              Vnew[i] = V[i] + dt * q[i] * Ex[i]/me;
          }
        #ifdef MODE_DEBUG
          for(i = 0; i < N; i++)
            {
            printf("Vnew = %.15lf \n",Vnew[i]);
            }
        #endif      
        // Updating Postion
        for(i = 0; i < N; i++)
          {
          Xnew[i] = X[i] + dt * Vnew[i];
          }
        #ifdef MODE_DEBUG  
          for(i = 0; i < N; i++)
            {
            printf("Xnew = %.15lf \n",Xnew[i]);
            } 
        #endif
        for(i = 0; i < N; i++)
          {
          X[i] = Xnew[i];
          V[i] = Vnew[i];
          }

    }
  for(i = 0; i < N; i++) {
    printf("Xnew = %.15lf \n",Xnew[i]);
  }
  // Free up the memory
  free(node);
  free(xi);
  free(xe);
  free(vi);
  free(ve);
  free(X);
  free(V);
  free(Vnew);
  free(Xnew);
 // free(hx);
  free(n);
//  free(w);
  free(q);
  free(qn);
  free(phi);

  // 'Point to data' must not be used again, unless re-assigned using malloc()
  node = NULL;
  xi   = NULL;
  xe   = NULL;
  vi   = NULL;
  ve   = NULL;
  X    = NULL;
  V    = NULL;
  Vnew = NULL;
  Xnew = NULL;
  //hx   = NULL;
  n    = NULL;
  //w    = NULL;
  q    = NULL;
  qn   = NULL;
  phi  = NULL;
  
  return 0;
}

double get_phi(double qn[])
{
  /* Locals */
  MKL_INT nn = NN, nrhs = NRHS, lda = LDA, ldb = LDB, info;
  /* Local arrays */
  MKL_INT ipiv[NN];
  double *aa   = (double *) malloc((NN*NN) * sizeof(double));
  if (aa == NULL) {
  printf("  Not enough memory to assign aa[%d]\n",NN*NN);
  exit(1);
  }else {
  printf("\n");
  printf("  Enough memory to assign aa[%d]\n",NN*NN);
  double aa[NN*NN];
  }
  int temp=1,i;
  for(i=0;i<NN*NN;i+=NN+temp)
  {
  aa[i]   = -2;
  aa[i+1] =  1;
  aa[i-1] =  1;
  }
  
  //double *bb   = (double *) malloc((NN) * sizeof(double));
  //if (bb == NULL) {
  //printf("  Not enough memory to assign bb[%d]\n",NN);
  //exit(1);
  //}else { 
  //printf("\n");
  //printf("  Enough memory to assign bb[%d]\n",NN);
  //double bb[NN];
  //}
  double e0 = 8.85;
  for(i=0;i<NN;i++)
    {
    bb[i] = 0.0;
    }
  bb[0] = qn[1]/e0;
  for(i=1;i<NN;i++)  
    {
    bb[i] = qn[i+1]/e0;
    }
  /* Executable statements */
  printf( "LAPACKE_dgesv (column-major, high-level) Example Program Results\n" );
  /* Solve the equations A*X = B */
  info = LAPACKE_dgesv( LAPACK_COL_MAJOR, nn, nrhs, aa, lda, ipiv,bb, ldb );
  /* Check for the exact singularity */
  if( info > 0 ) 
    {
    printf( "The diagonal element of the triangular factor of A,\n" );
    printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
    printf( "the solution could not be computed.\n" );
    exit( 1 );
    }
  #ifdef MODE_DEBUG  
  /* Print solution */
  print_matrix( "Solution", nn, nrhs, bb, ldb );
  /* Print details of LU factorization */
  print_matrix( "Details of LU factorization", nn, nn, aa, lda );
/* Print pivot indices */
  print_int_vector( "Pivot indices", nn, ipiv );
  #endif
  //exit( 0 );
  return 0;
} /* End of LAPACKE_dgesv Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT mm, MKL_INT nn, double* aa, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < mm; i++ ) {
                for( j = 0; j < nn; j++ ) printf( " %6.2lf", aa[i+j*lda] );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, MKL_INT nn, MKL_INT* aa ) {
        MKL_INT j;
        printf( "\n %s\n", desc );
        for( j = 0; j < nn; j++ ) printf( " %6i", aa[j] );
        printf( "\n" );
}


double get_field_pic(double phi[],double dn)
{
  int i;

  E[0]   = -( phi[1] - phi[0] )/dn;
  E[num_node-1] = -( phi[num_node-1] - phi[num_node-2] )/dn;
  
  for(i = 1; i < num_node-1; i++)
    {  
    E[i] = - ( phi[i+1] - phi[i-1] ) / ( 2 * dn );
    }
  return 0;
}

double interp_linear( double X[], double E[], double dn,double n[] )
{
    int i;
    double *slope   = (double *) malloc((num_node-1) * sizeof(double));
    if (slope == NULL) {
    printf("  Not enough memory to assign slope[%d]\n",num_node-1);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign slope[%d]\n",num_node-1);
    double slope[num_node-1];
    // memset( slope, 0.0, (num_node-1) * sizeof(double) );
    }
    for(i=0;i<num_node-1;i++)
      {
      slope[i]  = 0.0;
      }
    
    for(i = 0; i < num_node-1; i++)
      {
      slope[i] = (E[i+1]-E[i])/dn;
      }
    for(i = 0; i < N; i++)
      {
      Ex[i] = X[i]*slope[(int)(n[i]-1)] + E[(int)(n[i]-1)];
      }
  return 0;
}
