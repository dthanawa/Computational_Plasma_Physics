// Compile and Run the file with following command
// ./Direct_Summation.sh
// For Debug mode : gcc -Wall -g -fopenmp -DMODE_DEBUG Direct_Summation.c -lm -o Direct_Summation.x
// For Prod modde : gcc -Wall -g -fopenmp -DMODE_PROD Direct_Summation.c -lm -o Direct_Summation.x

#define PAR_OMP        // Conditional inclusion of mpi.h
#include<functions.h>
#include<stdio.h>
#include<math.h>

#define PI 3.141592653589793

#ifdef MODE_DEBUG
  #define Ni 10
  #define Ne 10
  #define c 0.01
#endif
#ifdef MODE_PROD
  #define Ni 100
  #define Ne 100
  #define c 0.01
#endif

int main(int argc, char **argv)
  {
  // Variable declaration

    int num_procs;           // Total number of available processors
    int max_threads;         // Maximum number of usable threads
    int thread_id;           // ID of each participating thread
    double wall_time = 0.00; // Wall time
  // Total number of available processors and maximum number of usable threads
    num_procs   = omp_get_num_procs();
    max_threads = omp_get_max_threads(); // Same as OMP_NUM_THREADS

  // Start the timer
    wall_time = omp_get_wtime();

  // MASTER thread
    thread_id = omp_get_thread_num();
    printf("\n");
    printf("  Total number of processors available : %d\n",   num_procs);
    printf("  Maximum number of usable threads     : %d\n",   max_threads);
    printf("  Thread ID of master processor        : %d\n\n", thread_id);
    
    int N = Ne + Ni;
    //int l=N;
  // Tempelate for memory preallocation
  /*
    double *a   = (double *) malloc((N) * sizeof(double));
    if (a == NULL) {
    printf("  Not enough memory to assign a[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign a[%d]\n", N);
    double a[N];
    memset( a, 0, N * sizeof(double) );
    }
  */
  // Memory Preallocation and Initialization BEGINS 


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
    double *yi   = (double *) malloc((N) * sizeof(double));
    if (yi == NULL) {
    printf("  Not enough memory to assign yi[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign yi[%d]\n", N);
    double yi[N];
    // memset( yi, 0.0, N * sizeof(double) );
    }
    double *Vix   = (double *) malloc((N) * sizeof(double));
    if (Vix == NULL) {
    printf("  Not enough memory to assign Vix[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Vix[%d]\n", N);
    double Vix[N];
    // memset( Vix, 0.0, N * sizeof(double) );
    }
    double *Viy   = (double *) malloc((N) * sizeof(double));
    if (Viy == NULL) {
    printf("  Not enough memory to assign Viy[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Viy[%d]\n", N);
    double Viy[N];
    // memset( Viy, 0.0, N * sizeof(double) );
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
    double *ye   = (double *) malloc((N) * sizeof(double));
    if (ye == NULL) {
    printf("  Not enough memory to assign ye[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign ye[%d]\n", N);
    double ye[N];
    // memset( ye, 0.0, N * sizeof(double) );
    }
    double *Vex   = (double *) malloc((N) * sizeof(double));
    if (Vex == NULL) {
    printf("  Not enough memory to assign Vex[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Vex[%d]\n", N);
    double Vex[N];
    // memset( Vex, 0.0, N * sizeof(double) );
    }
    double *Vey   = (double *) malloc((N) * sizeof(double));
    if (Vey == NULL) {
    printf("  Not enough memory to assign Vey[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Vey[%d]\n", N);
    double Vey[N];
    // memset( Vey, 0.0, N * sizeof(double) );
    }
    double *Y   = (double *) malloc((N) * sizeof(double));
    if (Y == NULL) {
    printf("  Not enough memory to assign Y[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Y[%d]\n", N);
    double Y[N];
    // memset( Y, 0.0, N * sizeof(double) );
    }
    double *Vx   = (double *) malloc((N) * sizeof(double));
    if (Vx == NULL) {
    printf("  Not enough memory to assign Vx[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Vx[%d]\n", N);
    double Vx[N];
    // memset( Vx, 0.0, N * sizeof(double) );
    }
    double *Vy   = (double *) malloc((N) * sizeof(double));
    if (Vy == NULL) {
    printf("  Not enough memory to assign Vy[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Vy[%d]\n", N);
    double Vy[N];
    // memset( Vy, 0.0, N * sizeof(double) );
    }    
    double *Ejx   = (double *) malloc((N) * sizeof(double));
    if (Ejx == NULL) {
    printf("  Not enough memory to assign Ejx[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Ejx[%d]\n", N);
    double Ejx[N];
    // memset( Ejx, 0.0, N * sizeof(double) );
    }
    double *Ejy   = (double *) malloc((N) * sizeof(double));
    if (Ejy == NULL) {
    printf("  Not enough memory to assign Ejy[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Ejy[%d]\n", N);
    double Ejy[N];
    // memset( Ejy, 0.0, N * sizeof(double) );
    }
    double *Ex   = (double *) malloc((N) * sizeof(double));
    if (Ex == NULL) {
    printf("  Not enough memory to assign Ex[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Ex[%d]\n", N);
    double Ex[N];
    // memset( Ex, 0.0, N * sizeof(double) );
    }    
    double *Ey   = (double *) malloc((N) * sizeof(double));
    if (Ey == NULL) {
    printf("  Not enough memory to assign Ey[%d]\n", N);
    exit(1);
    }else {
    printf("\n");
    printf("  Enough memory to assign Ey[%d]\n", N);
    double Ey[N];
    // memset( Ey, 0.0, N * sizeof(double) );
    }
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
    int i,j,k;

    for(i=0;i<N;i++)
      {
      xi[i]=0.0;
      yi[i]=0.0;
      xe[i]=0.0;
      ye[i]=0.0;
      X[i]=0.0;
      Y[i]=0.0;
      Vex[i]=0.0;
      Vey[i]=0.0;
      Ejx[i]=0.0;
      Ejy[i]=0.0;
      Ex[i]=0.0;
      Ey[i]=0.0;
      q[i]=0.0;
      Vx[i]=0.0;
      Vy[i]=0.0;
      
      }
     #ifdef MODE_DEBUG
    for(i=0;i<N;i++) {
      printf("xi[%d] = %lf \n",i,xi[i]);
    }
    #endif
    double d = 0.01, dt = 0.05;
    double dtheta, dgamma,t;
    double numx,numy,denxy;
 // Mass of ions and electrons
    double mi = 1000.0, me = 1.0;
        
    dtheta = 2.0*PI/Ni;
    dgamma = 2.0*PI/Ne;
 // Memory Preallocatiion and Initialization ENDS

    double temp = 1.0;

    for(i=0;i<Ni;i++)
        {
        xi[i]  = 1.0*cos((double)temp * dtheta);
        yi[i]  = 1.0*sin((double)temp * dtheta);
        Vix[i] = 0.0;
        Viy[i] = 0.0;
        temp++;
        }
    #ifdef MODE_DEBUG
      printf("i xi yi Vix Viy \n");  
      for(i=0;i<N;i++) {
        printf("%02d  % 017.15lf  % 017.15lf  % 017.15lf  % 017.15lf\n",i, xi[i],yi[i],Vix[i],Viy[i]);
      }
    #endif
    temp = 1.0;    
    for(i=Ni;i<N;i++)
        {
        xe[i]  = 1.02 * cos((double)temp * dgamma);
        ye[i]  = 1.02 * sin((double)temp * dgamma);
        Vex[i] = c    * cos((double)temp * dgamma);
        Vey[i] = c    * sin((double)temp * dgamma);
        temp++;
        }
    #ifdef MODE_DEBUG
      printf("xe ye Vex Vey\n");
      for(i=0;i<N;i++) {
        printf("%.15lf  %.15lf  %.15lf  %.15lf\n",xe[i],ye[i],Vex[i],Vey[i]);
      }
    #endif
    for(i=0;i<N;i++)
        {
        X[i]  =  xi[i] + xe[i];
        Y[i]  =  yi[i] + ye[i];
        Vx[i] = Vix[i] + Vex[i];
        Vy[i] = Viy[i] + Vey[i];
        }

    #ifdef MODE_DEBUG
      printf("X Y Vx Vy\n");
      for(i=0;i<N;i++) {
        printf("%.15lf  %.15lf  %.15lf  %.15lf\n",X[i],Y[i],Vx[i],Vy[i]);
      }
    #endif

    for(i=0;i<N;i++)
        {
        if(i<Ni)
        q[i]=1.0;
        else
        q[i]=-1.0;
        }

    #pragma omp parallel shared(max_threads) private(thread_id) 
        {
        // Thread ID
        thread_id = omp_get_thread_num();
        for(t=0; t<=20; t+=dt)
            {
            memset( Ex, 0, N * sizeof(double) );
            memset( Ey, 0, N * sizeof(double) );

            for(i=0;i<N;i++)
                {
                memset( Ejx, 0, N * sizeof(double) );
                memset( Ejy, 0, N * sizeof(double) );
                for(j=0;j<N;j++)
                    {
                    if(i!=j)
                      {
                      numx   = (1.0/(2*PI))*q[j]*(X[i]-X[j]);
                      numy   = (1.0/(2*PI))*q[j]*(Y[i]-Y[j]);
                      denxy  = sqrt(pow((X[i]-X[j]),2)+pow((Y[i]-Y[j]),2) + pow(d,2) );
                      Ejx[j] = numx/denxy;
                      Ejy[j] = numy/denxy;
                      }
                    else
                      {
                      Ejx[j] = 0;
                      Ejy[j] = 0;
                      }
                    }
                for(k=0;k<N;k++)
                    {
                    Ex[k]+=Ejx[k];
                    Ey[k]+=Ejy[k];
                    }
                }
        for(i=0;i<N;i++)
            {
            if(q[i]==1.0/Ni)
              {
              Vx[i] = Vx[i] + dt*q[i]*Ex[i]/mi;
              Vy[i] = Vy[i] + dt*q[i]*Ey[i]/mi;
              }
            else
              {
              Vx[i] = Vx[i] + dt*q[i]*Ex[i]/me;
              Vy[i] = Vy[i] + dt*q[i]*Ey[i]/me;
              }
            }
        for(i=0;i<N;i++)
            {
            X[i] = X[i] + dt*Vx[i];
            Y[i] = Y[i] + dt*Vy[i];
            }
   /*  
     for(i=0;i<N;i++)
        {
        printf("X at t(%.lf) =  %.15lf \n",t,X[i]);
        //printf("Y = %.lf \n",Y[i]); 
        }              
                
   */
            }
        }
  // Stop the timer and count the time
    wall_time = omp_get_wtime() - wall_time;

  // MASTER thread
    printf("\n");
    printf("  Total time taken : %.15lf seconds\n\n", wall_time);

    for(i=0;i<N;i++)
        {
        printf("X at t(%.lf) =  %.15lf \n",t,X[i]);
        //printf("Y = %.lf \n",Y[i]); 
        }

  // Free up the memory
  free(xi);
  free(yi);
  free(Vix);
  free(Viy);
  free(X);
  free(Y);
  free(xe);
  free(ye);
  free(Vex);
  free(Vey);
  free(Vx);
  free(Vy);
  free(Ejx);
  free(Ejy);
  free(Ex);
  free(Ey);
  free(q);
    
  // 'Point to data' must not be used again, unless re-assigned using malloc()
  xi  = NULL;  
  yi  = NULL;
  Vix = NULL; 
  Viy = NULL;  
  X   = NULL;  
  Y   = NULL;  
  xe  = NULL;  
  ye  = NULL;  
  Vex = NULL;
  Vey = NULL;
  Vx  = NULL;
  Vy  = NULL;
  Ejx = NULL;
  Ejy = NULL;
  Ex  = NULL;
  Ey  = NULL;
  q   = NULL;
    
    
    
  return 0;

}

