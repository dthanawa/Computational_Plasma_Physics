#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>

#define N 100
#define num_node 5
#define PI 3.141592653589793
//int N = 100;
//int num_node = 5;

// define the user-defined functions
// Electric Potential
double get_phi(double qn[]);

// Electric Field
double get_field_pic(double phi[],double dn);

// Linear Interpolation 
double interp_linear( double X[], double E[], double dn,double n[] );

// Global variables
//double *x;
//double *E;
//double *Ex;
//x = (double *)malloc((num_node-2) * sizeof(double));
//E = (double *)malloc(num_node * sizeof(double));
//Ex = (double *)malloc(N * sizeof(double));
//x[num_node-2]  = 0.0; E[num_node]  = 0.0; Ex[N] = 0.0;
//E = malloc(sizeof(*E) * num_node);
//Ex = malloc(sizeof(*Ex) * N);
//x = malloc(sizeof(*x) * (num_node-2));
double* E = 0; 
double* Ex = 0;
double* x = 0;
//x[num_node-2]  = 0.0; 
//E[num_node]  = 0.0;
//Ex[N] = 0.0;
// Main Function
int main()
{
    // Initialization of variables
    // #of ions and electrons

    int Ni = 50, Ne = 50;// N;
    int i, j;    
    double t, a1 = 0.0, b1 = 1.0;
    // double *node,*xi,*xe,*vi,*ve,*X,*V,*Vnew,*Xnew,**hx,*n,**w,*q,*qn,*phi;
   
    // Memory pre-alocation and iff successful, initialization
// Template
// If there isn't enough memory, display an appropriate message and exit
// If there is enough memory, declare and initialize a ZERO array
// If need be, display the array elements and confirm that they are ZEROs
// double *x_ptr = (double *) malloc(N * sizeof(double));
// if (x_ptr == NULL) {
//   printf("\n");
//   printf("  Not enough memory to assign x[%d]\n", N);
//   exit(1);
// } else {
//   printf("\n");
//   printf("  Enough memory to assign x[%d]\n", N);
//   double x[N] = {0};
//
//   // DEBUG
//   // Confirm that all elements of x[N] are 0.0
//   printf("\n");
//   for (i = 0; i < N; i++) {
//     printf("    x[%d] = %lf\n", i, x[i]);
//   }
// }

    printf("\n");
    double *ptr_node = (double *) malloc(num_node * sizeof(double));
    if (ptr_node == NULL) {
      printf("  Not enough memory to assign node[%d]\n", num_node);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign node[%d]\n", num_node);
      double node[num_node] = {0};

      // DEBUG
      // Confirm that all elements of node[num_node] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < num_node; i++) {
          printf("    node[%d] = %lf\n", i, node[i]);
        }
      #endif
    }

    double *ptr_xi   = (double *) malloc(N * sizeof(double));
    if (ptr_xi == NULL) {
      printf("  Not enough memory to assign xi[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign xi[%d]\n", N);
      double xi[N] = {0};

      // DEBUG
      // Confirm that all elements of xi[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    xi[%d] = %lf\n", i, xi[i]);
        }
      #endif
    }
    printf("\n");
    double *ptr_xe = (double *) malloc(N * sizeof(double));
    if (ptr_xe == NULL) {
      printf("  Not enough memory to assign xe[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign xe[%d]\n", N);
      double xe[N] = {0};

      // DEBUG
      // Confirm that all elements of xe[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    xe[%d] = %lf\n", i, xe[i]);
        }
      #endif
    }
    double *ptr_vi   = (double *) malloc(N * sizeof(double));
    if (ptr_vi == NULL) {
      printf("  Not enough memory to assign vi[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign vi[%d]\n", N);
      double vi[N] = {0};

      // DEBUG
      // Confirm that all elements of xi[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    vi[%d] = %lf\n", i, vi[i]);
        }
      #endif
    }
    double *ptr_ve   = (double *) malloc(N * sizeof(double));
    if (ptr_ve == NULL) {
      printf("  Not enough memory to assign ve[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign ve[%d]\n", N);
      double ve[N] = {0};

      // DEBUG
      // Confirm that all elements of vei[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    ve[%d] = %lf\n", i, ve[i]);
        }
      #endif
    }
    double *ptr_X   = (double *) malloc(N * sizeof(double));
    if (ptr_X == NULL) {
      printf("  Not enough memory to assign X[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign X[%d]\n", N);
      double X[N] = {0};

      // DEBUG
      // Confirm that all elements of X[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    X[%d] = %lf\n", i, X[i]);
        }
      #endif
    }
    double *ptr_V   = (double *) malloc(N * sizeof(double));
    if (ptr_V == NULL) {
      printf("  Not enough memory to assign V[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign V[%d]\n", N);
      double V[N] = {0};

      // DEBUG
      // Confirm that all elements of V[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    V[%d] = %lf\n", i, V[i]);
        }
      #endif
    }
    double *ptr_Vnew   = (double *) malloc(N * sizeof(double));
    if (ptr_Vnew == NULL) {
      printf("  Not enough memory to assign Vnew[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign Vnew[%d]\n", N);
      double Vnew[N] = {0};

      // DEBUG
      // Confirm that all elements of Vnew[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    Vnew[%d] = %lf\n", i, Vnew[i]);
        }
      #endif
    }
    double *ptr_Xnew   = (double *) malloc(N * sizeof(double));
    if (ptr_Xnew == NULL) {
      printf("  Not enough memory to assign Xnew[%d]\n", N);
      exit(1);
    } else {
      printf("\n");  
      printf("  Enough memory to assign Xnew[%d]\n", N);
      double Xnew[N] = {0};

      // DEBUG
      // Confirm that all elements of Xnew[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    Vnew[%d] = %lf\n", i, Xnew[i]);
        }
      #endif
    } 
   double *ptr_hx   = (double *) malloc(N * sizeof(double));
    if (ptr_hx == NULL) {
      printf("  Not enough memory to assign hx[%d]\n", N);
      exit(1);
    } else {
      printf("\n");  
      printf("  Enough memory to assign hx[%d]\n", N);
      double hx[N][2] = {{0}};
/*
      // DEBUG
      // Confirm that all elements of Vnew[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    hx[%d] = %lf\n", i, hx[i]);
        }
      #endif
*/
    } 
    double *ptr_n   = (double *) malloc(N * sizeof(double));
    if (ptr_n == NULL) {
      printf("  Not enough memory to assign n[%d]\n", N);
      exit(1);
    } else {
      printf("\n");  
      printf("  Enough memory to assign n[%d]\n", N);
      double n[N] = {0};

      // DEBUG
      // Confirm that all elements of n[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    n[%d] = %lf\n", i, n[i]);
        }
      #endif
    }
   double *ptr_w   = (double *) malloc(N * sizeof(double));
    if (ptr_w == NULL) {
      printf("  Not enough memory to assign w[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign w[%d]\n", N);
      double w[N][2] = {{0}};
/*
      // DEBUG
      // Confirm that all elements of Vnew[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    hx[%d] = %lf\n", i, hx[i]);
        }
      #endif
*/
    }
 
    double *ptr_q   = (double *) malloc(N * sizeof(double));
    if (ptr_q == NULL) {
      printf("  Not enough memory to assign q[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign q[%d]\n", N);
      double q[N] = {0};

      // DEBUG
      // Confirm that all elements of n[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    q[%d] = %lf\n", i, q[i]);
        }
      #endif
    }

    

    
    double *ptr_qn = (double *) malloc(num_node * sizeof(double));
    if (ptr_qn == NULL) {
      printf("  Not enough memory to assign qn[%d]\n", num_node);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign qn[%d]\n", num_node);
      double qn[num_node] = {0};

      // DEBUG
      // Confirm that all elements of node[num_node] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < num_node; i++) {
          printf("    qn[%d] = %lf\n", i, qn[i]);
        }
      #endif
    }

    double *ptr_phi   = (double *) malloc(N * sizeof(double));
    if (ptr_phi == NULL) {
      printf("  Not enough memory to assign n[%d]\n", N);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign n[%d]\n", N);
      double phi[N] = {0};

      // DEBUG
      // Confirm that all elements of n[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    phi[%d] = %lf\n", i,phi[i]);
        }
      #endif
    }


    printf("\n");

    double dxe, dxi, m = 1.0, k = 0.0, mi = 1000.0, me = 1.0;
    
    // Initialization
   
     double xi[N]          = {0};
    double xe[N]          = {0};
    double vi[N]          = {0};
    double ve[N]          = {0};
    double X[N]           = {0};
    double V[N]           = {0};
    double Vnew[N]        = {0};
    double Xnew[N]        = {0};
    double hx[N][2]       = {{0}};
    double n[N]           = {0};
    double w[N][2]        = {{0}};
    double q[N]           = {0};
    double qn[num_node]   = {0};
    double phi[num_node]  = {0};
    double node[num_node] = {0};
   
    double dt             = 0.0005;
    
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
    
    // Merege Velocity and Postion of ions and electrons in one vector    
    for(i = 0; i < N; i++)
        {
            X[i] = xi[i] + xe[i];
            V[i] = vi[i] + ve[i];
            printf("X = %.15lf\t\t V = %.15lf\n",X[i],V[i]);
        }
    
    // Creating 1D Mesh   
    for(i = 0; i < num_node; i++)
        {
            node[i] = k;
            printf("node(%d) = %.15lf\n",i,node[i]);
            k = k + 1.0/(num_node-1.0);
        }
        
    double dn = 1.0/(num_node-1.0);
    
    // Time Loop BEGINS
    for(t = 0; t <= 1.5; t+=dt)
        {
            // Calculating hx and Weight fraction
            for(i = 0;i<N;i++)
                {
                    //n[i] = ceil(X[i]/dn);
                    n[i] = floor((X[i] - a1)/(b1-a1)*(num_node - 1)+1);
                    hx[i][1] =  X[i] - (n[i]-1.0)*dn;
                    hx[i][2] = n[i]*dn - X[i];
                    printf("n = %.15lf hx = %.15lf  %.15lf\n",n[i],hx[i][1],hx[i][2]);
                }
        
            for(i = 1;i<N;i++)
                {
                    w[i][1] = hx[i][1]/(hx[i][1]+hx[i][2]);
                    w[i][2] = 1.0-w[i][1];
                    printf("w1 = %.15lf  w2 = %.15lf\n",w[i][1],w[i][2]);
                }
            
            // Random Charge generator
            srand((unsigned int)time(NULL));
            
            for(i=0;i<N/2;i++)
                q[i] = ((double)rand()/(double)(RAND_MAX)) * r;
        
            for(i=N/2;i<N;i++)
                q[i] = ((double)rand()/(double)(RAND_MAX)) * (-r);
    
            for(i = 1;i<N;i++)
                printf("q = %.15lf \n",q[i]);
            
            // Calculating Charge density
            for(j = 1;j<N;j++)
                {
                    n[j] = ceil(X[j]/dn);
                    qn[(int)n[j]] = qn[(int)n[j]] +  w[j][2]*q[j];
                    qn[(int)n[j]+1] = qn[(int)n[j]+1] +  w[j][1]*q[j];
                }
            for(i=0;i<num_node;i++)
                printf("qn = %.15lf \n",qn[i]);
            
            // Calling Electric Potential Function
            get_phi(qn);
            
            for(i=1;i<num_node-1;i++)
                phi[i] = x[i-1];
    
            // Calling Electric Field Function
            get_field_pic(phi,dn);
        
            for(i = 0;i < num_node; i++)
                printf("phi = %.15lf\t E = %.15lf \n",phi[i],E[i]);
            
            // Calling Linear Interpolation Function
            interp_linear( X, E, dn, n );
        
            for(i = 0; i < N; i++)
                printf("Ex = %.15lf \n",Ex[i]);
            
            // Updating Velocity
            for(i = 0;i < N; i++)
                {
                    if(q[i] > 0)
                        Vnew[i] = V[i] + dt * q[i] * Ex[i]/mi;
                    else
                        Vnew[i] = V[i] + dt * q[i] * Ex[i]/me;
                }
            
            for(i = 0; i < N; i++)
                printf("Vnew = %.15lf \n",Vnew[i]);
            
            // Updating Postion
            for(i = 0; i < N; i++)
                Xnew[i] = X[i] + dt * Vnew[i]; 
    
            for(i = 0; i < N; i++)
                printf("Xnew = %.15lf \n",Xnew[i]);
    
            for(i = 0; i < N; i++)
                {
                    X[i] = Xnew[i];
                    V[i] = Vnew[i];
                }
            
            printf("t = %.lf\n",t);
            
        }    
    

  // Free up the memory
  free(ptr_node);
  free(ptr_xi);
  free(ptr_xe);
  free(ptr_vi);
  free(ptr_ve);
  free(ptr_X);
  free(ptr_V);
  free(ptr_Vnew);
  free(ptr_Xnew);
  free(ptr_hx);
  free(ptr_n);
  free(ptr_w);
  free(ptr_q);
  free(ptr_qn);
  free(ptr_phi);

  // 'Point to data' must not be used again, unless re-assigned using malloc()
  ptr_node = NULL;
  ptr_xi   = NULL;
  ptr_xe   = NULL;
  ptr_vi   = NULL;
  ptr_ve   = NULL;
  ptr_X    = NULL;
  ptr_V    = NULL;
  ptr_Vnew = NULL;
  ptr_Xnew = NULL;
  ptr_hx   = NULL;
  ptr_n    = NULL;
  ptr_w    = NULL;
  ptr_q    = NULL;
  ptr_qn   = NULL;
  ptr_phi  = NULL;

    return 0;
}

// Electric Potential Function
double get_phi(double qn[])
{
    int tn, i, j;
    double e0 = 8.85;
    tn = num_node-2;
   // double **a,*f,*alfa,*bet;
        
    double *ptr_a   = (double *) malloc((num_node-2) * sizeof(double));
    if (ptr_a == NULL) {
      printf("  Not enough memory to assign a[%d]\n", num_node-2);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign a[%d]\n", num_node-2);
      double a[num_node-2][num_node-2] = {{0}};
/*
      // DEBUG
      // Confirm that all elements of Vnew[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < N; i++) {
          printf("    hx[%d] = %lf\n", i, hx[i]);
        }
      #endif
*/
    }
   // a = (double **)malloc((num_node-2) * sizeof(double));
   
     double *ptr_f   = (double *) malloc((num_node-2) * sizeof(double));
    if (ptr_f == NULL) {
      printf("  Not enough memory to assign f[%d]\n", num_node-2);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign f[%d]\n", num_node-2);
      double f[num_node-2] = {0};

      // DEBUG
      // Confirm that all elements of n[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < num_node-2; i++) {
          printf("    f[%d] = %lf\n", i,f[i]);
        }
      #endif
    }

   // f = (double *)malloc((num_node-2) * sizeof(double));
    
    double *ptr_alfa   = (double *) malloc((num_node-2) * sizeof(double));
    if (ptr_alfa == NULL) { 
      printf("  Not enough memory to assign alfa[%d]\n", num_node-2);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign alfa[%d]\n", num_node-2);
      double alfa[num_node-2] = {0};

      // DEBUG
      // Confirm that all elements of n[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < num_node-2; i++) {
          printf("    alfa[%d] = %lf\n", i,alfa[i]);
        }
      #endif
    }
    
double *ptr_bet   = (double *) malloc((num_node-2) * sizeof(double));
    if (ptr_bet == NULL) { 
      printf("  Not enough memory to assign bet[%d]\n", num_node-2);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign bet[%d]\n", num_node-2);
      double bet[num_node-2] = {0};

      // DEBUG
      // Confirm that all elements of n[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < num_node-2; i++) {
          printf("    bet[%d] = %lf\n", i,bet[i]);
        }
      #endif
    }   
    //bet = (double *)malloc((num_node-2) * sizeof(double));
    //a[num_node-2][num_node-2] = 0.0; f[num_node-2] = 0.0; alfa[num_node-2] = 0.0; bet[num_node-2] = 0.0;
  
  //Initialization
  double a[num_node-2][num_node-2] = {{0}};
   double f[num_node-2] = {0}; 
   double alfa[num_node-2] = {0}; 
   double bet[num_node-2] = {0}; 
    
    
    a[0][0] = -2;
    a[0][1] = 1;
    f[0] = qn[1]/e0;
    
    for(i=1;i<tn;i++)
        {
            a[i][i-1] = 1;
            a[i][i]   = -2;
            if(i!=num_node-3)
            a[i][i+1] = 1;
            f[i] = qn[i+1]/e0; 
        }
        
    for(i=0;i<tn;i++)
        {
           for(j=0;j<tn;j++)
                printf("%.15lf ",a[i][j]);
        
            printf("\n");
        }
        
    alfa[0] = a[0][0];
    bet[0]  = f[0]/a[0][0];
    
    for(i=1;i<num_node-2;i++)
    {
        alfa[i] = a[i][i] - (a[i][i-1])*a[i-1][i]/alfa[i-1];
        bet[i] = (f[i] - (a[i][i-1])*bet[i-1])/alfa[i];
    }
    
    x[num_node-3] = bet[num_node-3];
    
    for(i = num_node-4; i >= 0; i--)
        x[i] = bet[i]-(a[i][i+1])*x[i+1]/alfa[i];
    
    for(i = 0;i<3;i++)
        printf("x=%f\n",x[i]);
    
    return 0;
}
double get_field_pic(double phi[],double dn)
{
    int i;

    E[0]   = -( phi[1] - phi[0] )/dn;
    E[num_node-1] = -( phi[num_node-1] - phi[num_node-2] )/dn;

    for(i = 1; i < num_node-1; i++)
        E[i] = - ( phi[i+1] - phi[i-1] ) / ( 2 * dn );

    return 0;
}

double interp_linear( double X[], double E[], double dn,double n[] )
{
    int i;
    //double *slope;
    int tn = num_node-1;
    double *ptr_slope   = (double *) malloc((num_node-1) * sizeof(double));
    if (ptr_slope == NULL) { 
      printf("  Not enough memory to assign slope[%d]\n", num_node-1);
      exit(1);
    } else {
      printf("\n");
      printf("  Enough memory to assign slope[%d]\n", num_node-1);
      double slope[num_node-1] = {0};

      // DEBUG
      // Confirm that all elements of n[N] are 0.0
      #ifdef MODE_DEBUG
        printf("\n");
        for (i = 0; i < num_node-1; i++) {
          printf("    slope[%d] = %lf\n", i,slope[i]);
        }
      #endif
    }
  //  slope = (double *)malloc((num_node-1) * sizeof(double));
    // slope[num_node-1] = 0.0;
    double slope[num_node-1] = {0}; 
    for(i = 0; i < num_node-1; i++)
        slope[i] = (E[i+1]-E[i])/dn;
    
    for(i = 0; i < N; i++)
        Ex[i] = X[i]*slope[(int)(n[i]-1)] + E[(int)(n[i]-1)];
    return 0;
}
