// Headers
#define PAR_OMP        // Conditional inclusion of mpi.h
#include <functions.h>

#include <stdio.h>
#include<math.h>
#ifdef MODE_DEBUG
  #define l 2
#endif
#ifdef MODE_PROD
  #define l 2000
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
  printf("\n");
  printf("  Total number of processors available : %d\n",   num_procs);
  printf("  Maximum number of usable threads     : %d\n\n", max_threads);

        int Ni = 1000, Ne = 1000,i,j,k;
        int N = Ne + Ni;
    //const int l = 108;
        float xi[l] = { 0.0 },yi[l] = { 0.0 },Vix[l] = { 0.0 },Viy[l] = { 0.0 },X[l] = { 0.0 };
        float xe[l] = { 0.0 },ye[l] = { 0.0 },Vex[l] = { 0.0 },Vey[l] = { 0.0 },Y[l] = { 0.0 };
        float dtheta, dgamma , c=0.01, mi=1000.0, me = 1.0,t;
        float Vx[l] = { 0.0 }, Vy[l] = { 0.0 }, Ejx[l] = { 0.0 },Ejy[l] = { 0.0 };
        float Ex[l] = { 0.0 }, Ey[l] = { 0.0 }, q[l] = {0.0}, d = 0.01, dt = 0.05;
        dtheta = 2.0*(22.0/7.0)/Ni;
        dgamma = 2.0*(22.0/7.0)/Ne;
        for(i=0;i<Ni;i++)
            {
                xi[i] = cos(i*dtheta);
                yi[i] = sin(i*dtheta);
                Vix[i] = 0;
                Viy[i] = 0;
            }
        for(i=8;i<Ne+8;i++)
            {
                xe[i] = 1.02*cos(i*dgamma);
                ye[i] = 1.02*sin(i*dgamma);
                Vex[i] = c*cos(i*dgamma);
                Vey[i] = c*sin(i*dgamma);
            }
        
        for(i=0;i<N;i++)
            {
                X[i]  = xi[i]+xe[i];
                Y[i]  = yi[i]+ye[i];
                Vx[i] = Vix[i]+Vex[i];
                Vy[i] = Viy[i]+Vey[i];
            }
#ifdef MODE_DEBUG
        for(i=0;i<N;i++)
            {
                // printf("%f\t%f\t%f\t%f\n",X[i],Y[i],Vx[i],Vy[i]);
                printf("%10.8f  %10.8f  %10.8f  %10.8f\n",X[i],Y[i],Vx[i],Vy[i]);
            }
#endif
        for(i=0;i<N;i++)
            {
                if(i<Ni)
                q[i]=1.0;
                else
                q[i]=-1.0;
            }
        for(t=0; t<=20; t+=dt)
            {
                float Ex[l] = { 0.0 },Ey[l] = { 0.0 };
#pragma omp parallel shared(max_threads) private(thread_id) 
  {
    // Thread ID
    thread_id = omp_get_thread_num();

    // Print a message
//    printf("    Hello, World! from thread %d out of %d\n", thread_id, max_threads);
  

                for(i=0;i<N;i++)
                    {
                        float Ejx[l] = { 0.0 },Ejy[l] = { 0.0 };
                        for(j=0;j<N;j++)
                            {
                                if(i!=j)
                                    {
                                        float numx,numy,denxy;
                                        numx   = (1/(2*(22/7)))*q[j]*(X[i]-X[j]);
                                        numy   = (1/(2*(22/7)))*q[j]*(Y[i]-Y[j]);
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
      for(i=0;i<N;i++)
                    {
                       printf("X at t(%f) =  %10.8f \n",t,X[i]);
                      //printf("Y = %.f \n",Y[i]); 
                    }              
                
   }        
        
  // Stop the timer and count the time
  wall_time = omp_get_wtime() - wall_time;

  // MASTER thread
  printf("\n");
  printf("  Total time taken : %10.8f seconds\n\n", wall_time);      
        
        
        return 0;   
        
    }
