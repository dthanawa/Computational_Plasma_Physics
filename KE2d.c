#include <stdio.h>
#include<math.h>
int main()
    {
        int Ni = 8, Ne = 100,N,i,j,k;
        N = Ne + Ni;
        float xi[108] = { 0.0 },yi[108] = { 0.0 },Vix[108] = { 0.0 },Viy[108] = { 0.0 },X[108] = { 0.0 };
        float xe[108] = { 0.0 },ye[108] = { 0.0 },Vex[108] = { 0.0 },Vey[108] = { 0.0 },Y[108] = { 0.0 };
        float dtheta, dgamma , c=0.01, mi=1000.0, me = 1.0,t;
        float Vx[108] = { 0.0 }, Vy[108] = { 0.0 }, Ejx[108] = { 0.0 },Ejy[108] = { 0.0 };
        float Ex[108] = { 0.0 }, Ey[108] = { 0.0 }, q[108] = {0.0}, d = 0.01, dt = 0.05;
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
        for(i=0;i<N;i++)
            {
                printf("%f\t%f\t%f\t%f\n",X[i],Y[i],Vx[i],Vy[i]);
            }
        for(i=0;i<N;i++)
            {
                if(i<Ni)
                q[i]=1.0;
                else
                q[i]=-1.0;
            }
        for(t=0; t<=20; t+=dt)
            {
                float Ex[108] = { 0.0 },Ey[108] = { 0.0 };
                for(i=0;i<N;i++)
                    {
                        float Ejx[108] = { 0.0 },Ejy[108] = { 0.0 };
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
                       printf("X at t(%f) =  %f \n",t,X[i]);
                      //printf("Y = %.f \n",Y[i]); 
                    }              
                
            }        
        
        
        
        
        return 0;   
        
    }
