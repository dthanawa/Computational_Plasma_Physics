#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#define PI (3.141592653589793)
float * get_phi((float) qn[5],(int) num_node);

int main()
{
     int Ni = 8, Ne = 8,i,j,N,num_node =5;
     N = Ni + Ne;
     float dxe,dxi,m=1.0,node[5]={ 0.0 },k=0.0;
     float xi[16] = { 0.0 },xe[16] = { 0.0 },vi[16] = { 0.0 },ve[16] = { 0.0 },X[16] = { 0.0 },V[16] = { 0.0 };
     float hx[16][2] = { 0.0 },n[16]={ 0.0 },w[16][2] = { 0.0 },q[16] = { 0.0 },qn[5] = { 0.0 },*phi;
     srand((unsigned int)time(NULL));
     float r = 2;
     dxe = 2.0/Ne;
     dxi = 1.0/Ni;
     
     for(i = 0 ; i< Ni;i++)
        {
            xi[i] = (m-0.5)*dxi;
            vi[i] = 0;
            m++;
        }
    m=1.0;
    for(i=Ne;i<3*Ne/2;i++)
        {
            xe[i] = (m-0.5)*dxe;
            ve[i] = 0.5 + 0.1 * sin(2.0*PI*xe[i]);
            m++;
        }
    for(i=3*Ne/2;i<2*Ne;i++)
        {
            xe[i] = xe[i-Ne/2];
            ve[i] = -0.5 - 0.1*sin(2.0*PI*xe[i]);
        }        
    for(i=0;i<N;i++)
        {
            X[i]=xi[i]+xe[i];
            V[i]=vi[i]+ve[i];
            printf("X = %.5f\t\t V = %.5f\n",X[i],V[i]);
        }
       
    for(i=0;i<num_node;i++)
        {
            node[i] = k;
            printf("node(%d) = %.5f\n",i,node[i]);
            k= k + 1.0/(num_node-1.0);
        }
    float dn =  1.0/(num_node-1.0);
    
    for(i = 0;i<N;i++)
        {
            n[i] = ceil(X[i]/dn);
            hx[i][1] =  X[i] - (n[i]-1.0)*dn;
            hx[i][2] = n[i]*dn - X[i];
            printf("hx = %.5f  %.5f\n",hx[i][1],hx[i][2]);
        }
    for(i = 1;i<N;i++)
        {
            w[i][1] = hx[i][1]/(hx[i][1]+hx[i][2]);
            w[i][2] = 1.0-w[i][1];
            printf("w1 = %.5f  w2 = %.5f\n",w[i][1],w[i][2]);
        }
        for(i=0;i<N/2;i++)
        {
            q[i] = ((float)rand()/(float)(RAND_MAX)) * r;
        }
        for(i=N/2;i<N;i++)
        {
            q[i] = ((float)rand()/(float)(RAND_MAX)) * (-r);
        }
        for(i = 1;i<N;i++)
        {
            printf("q = %.5f \n",q[i]);
        }
        for(j = 1;j<N;j++)
        {
            n[j] = ceil(X[j]/dn);
            qn[(int)n[j]] = qn[(int)n[j]] +  w[j][2]*q[j];
            qn[(int)n[j]+1] = qn[(int)n[j]+1] +  w[j][1]*q[j];
        }
        for(i=0;i<num_node;i++)
        {
            printf("qn = %.5f \n",qn[i]);
        }
        phi = get_phi((float) qn[5],(int) num_node);
    return 0;
}

float * get_phi((float) qn[5],(int) num_node)
{
    int n,i;
    float e0 = 8.85;
    n = num_node-2;
    float a[3][3] = { 0.0 }, f[3] = { 0.0 }, alfa[3] = { 0.0 }, bet[3] = { 0.0 },x[3] = { 0.0 };
    a[0][0] = -2;
    a[0][1] = 1;
    f[0]=qn[1]/e0;
    
    for(i=1;i<n;i++)
    {
        a[i][i-1]=1;
        a[i][i]=-2;
        a[i][i+1]=1;
        f[i]=qn[i+1]/e0; 
    }
    alfa[0]=a[0][0];
    bet[0]=f[0]/a[1][1];
    for(i=1;i<n;i++)
    {
        alfa[i] = a[i][i]-a[i][i-1]*a[i-1][i]/alfa[i-1];
        bet[i] = (f[i]-a[i][i-1]*bet[i-1])/alfa[i];
    }
    x[n-1]=bet[n-1];
    for(i=n-2;i>=0;i--)
    {
        x[i]=bet[i]-a[i][i+1]*x[i+1]/alfa[i];
    }
    
    
    return x;

}


