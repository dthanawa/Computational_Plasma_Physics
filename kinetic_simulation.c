#include<stdio.h>
#include<math.h>

int main()
    {
    double dxi,dxe;
    double xe[1000],xi[1000],ve[1000],vi[1000];
    double X[1000],V[1000],E[1000],Ej[1000];
    int Ne=4,Ni=4,N,i,t,k,j;
    double m=1.0;
    signed int qi[1000],qe[1000],q[1000];
    N=Ne+Ni;
    dxe=2.0/Ne;
    dxi=1.0/Ni;
    for(i=0;i<Ni;i++)
        {
        xi[i] = ((m-0.5)*dxi);
        vi[i] = 0;
        m++;
        printf("xi = %f  \n",xi[i]);
        }
    printf("x(1)=%f\n\n",xi[1]);
    m=1.0;
    for(i=0;i<Ne/2;i++)
        {
        xe[i] = (m-0.5)*dxe;
        ve[i] = 0.5 + 0.1*sin(2*22/7*xe[i]);
        m++;
        printf("xe(%d) = %f  \n",i,xe[i]);
        }
    for(i=(Ne/2);i<Ne;i++)
        {
        xe[i] = xe[i-(Ne/2)];
        ve[i] = -0.5 - 0.1*sin(2*22/7*xe[i]);
        printf("xe(%d) = %f  \n",i,xe[i]);
        }
    for(i=0;i<Ni;i++)
        {
        qi[i]=1;
        }
    for(i=0;i<Ne;i++)
        {
        qe[i]=-1;
        }
    for(i=0;i<Ni;i++)
        {
        V[i]=vi[i];
        X[i]=xi[i];
        q[i]=qi[i];
        //printf("v = %d \n\n V= %d \n",vi[i],V[i]);
        }
    for(i=Ni;i<N;i++)
        {
        V[i]=ve[i-Ne];
        X[i]=xe[i-Ne];
        q[i]=qe[i-Ne];
        }    
    for(i=0;i<N;i++)
        {
        printf("V= %f \t\t X= %f \t\t q= %d \n",V[i],X[i],q[i]);
        }
    int mi = 1000, me = 1,i1,i2;
    double dt=5*pow(10,-2),Vnew[1000],Xnew[1000],sum=0;
    //for(t=0;t<1.5;t+=dt)
       // {
        for(i=0;i<N;i++)
            {
            for(i1=0;i1<N;i1++)    
                {
                Ej[i1]=0;
                }
            for(j=0;j<N;j++)
                {
                if(X[i]>X[i])
                    Ej[j]=(q[j])/(-2.0);
                else if(X[i]<X[j])
                    Ej[j]=(q[j])/(2.0);
                else Ej[j]=0;
                }
                sum=0;
                for(i2=0;i2<N;i2++)
                {
                    sum+=Ej[i2];
                }
                E[i]=sum;
            }
                    
        for(i=0;i<N;i++)
            {
            if(q[i]==1)
                {
                Vnew[i] = V[i] + dt*(q[i])*(E[i])/mi;
                }
            else
                {
                Vnew[i] = V[i] + dt*(q[i])*(E[i])/me;
                }
                
            }
        for(i=0;i<N;i++)
            {
            Xnew[i] = X[i] + dt*V[i];
            printf("Xnew(%d) = %f  \n",i,Xnew[i]);
            }
            
      //  }
    }




