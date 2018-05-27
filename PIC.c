#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#define PI (3.141592653589793)
#define N (100)
#define num_node (5)
// define the user-defined functions
// Electric Potential
float get_phi(float qn[]);

// Electric Field
float get_field_pic(float phi[],float dn);

// Linear Interpolation 
float interp_linear( float X[], float E[], float dn,float n[] );

// Global variables
float x[num_node-2]  = {{ 0.0} }, E[num_node]  = { {0.0} }, Ex[N] = { {0.0} };

// Main Function
int main()
{
    // Initialization of variables
    // #of ions and electrons
    int Ni = 50, Ne = 50;// N;
    //N = Ni + Ne;

    // #of nodes
    //int num_node = 5; //p = 5;

    int i, j;    
    float t, a = 0.0, b = 1.0;
    float dxe, dxi, m = 1.0, node['p'] = { 0.0 }, k = 0.0, mi = 1000.0, me = 1.0;
    
    // Memory Preallocation for Velocity, postion, weight fraction, Charge density and Electric Potential
    float xi[N] = { 0.0 }, xe[N] = { 0.0 }, vi[N] = { 0.0 }, ve[N] = { 0.0 };
    float X[N] = { 0.0 }, V[N] = { 0.0 }, Vnew[N] = { 0.0 }, Xnew[N] = { 0.0 };
    float hx[N][2] = { 0.0 }, n[N] = { 0.0 }, w[N][2] = { 0.0 };
    float q[N] = { 0.0 },qn[num_node] = { 0.0 }, phi[num_node] = { 0.0 };
    
    // Time step
    float dt = 0.0005;
    
    float r = 2;
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
            printf("X = %.5f\t\t V = %.5f\n",X[i],V[i]);
        }
    
    // Creating 1D Mesh   
    for(i = 0; i < num_node; i++)
        {
            node[i] = k;
            printf("node(%d) = %.5f\n",i,node[i]);
            k = k + 1.0/(num_node-1.0);
        }
        
    float dn = 1.0/(num_node-1.0);
    
    // Time Loop BEGINS
    for(t = 0; t <= 1.5; t+=dt)
        {
            // Calculating hx and Weight fraction
            for(i = 0;i<N;i++)
                {
                    //n[i] = ceil(X[i]/dn);
                    n[i] = floor((X[i] - a)/(b-a)*(num_node - 1)+1);
                    hx[i][1] =  X[i] - (n[i]-1.0)*dn;
                    hx[i][2] = n[i]*dn - X[i];
                    printf("n = %.5f hx = %.5f  %.5f\n",n[i],hx[i][1],hx[i][2]);
                }
        
            for(i = 1;i<N;i++)
                {
                    w[i][1] = hx[i][1]/(hx[i][1]+hx[i][2]);
                    w[i][2] = 1.0-w[i][1];
                    printf("w1 = %.5f  w2 = %.5f\n",w[i][1],w[i][2]);
                }
            
            // Random Charge generator
            srand((unsigned int)time(NULL));
            
            for(i=0;i<N/2;i++)
                q[i] = ((float)rand()/(float)(RAND_MAX)) * r;
        
            for(i=N/2;i<N;i++)
                q[i] = ((float)rand()/(float)(RAND_MAX)) * (-r);
    
            for(i = 1;i<N;i++)
                printf("q = %.5f \n",q[i]);
            
            // Calculating Charge density
            for(j = 1;j<N;j++)
                {
                    n[j] = ceil(X[j]/dn);
                    qn[(int)n[j]] = qn[(int)n[j]] +  w[j][2]*q[j];
                    qn[(int)n[j]+1] = qn[(int)n[j]+1] +  w[j][1]*q[j];
                }
            for(i=0;i<num_node;i++)
                printf("qn = %.5f \n",qn[i]);
            
            // Calling Electric Potential Function
            get_phi(qn);
            
            for(i=1;i<num_node-1;i++)
                phi[i] = x[i-1];
    
            // Calling Electric Field Function
            get_field_pic(phi,dn);
        
            for(i = 0;i < num_node; i++)
                printf("phi = %.5f\t E = %.5f \n",phi[i],E[i]);
            
            // Calling Linear Interpolation Function
            interp_linear( X, E, dn, n );
        
            for(i = 0; i < N; i++)
                printf("Ex = %.5f \n",Ex[i]);
            
            // Updating Velocity
            for(i = 0;i < N; i++)
                {
                    if(q[i] > 0)
                        Vnew[i] = V[i] + dt * q[i] * Ex[i]/mi;
                    else
                        Vnew[i] = V[i] + dt * q[i] * Ex[i]/me;
                }
            
            for(i = 0; i < N; i++)
                printf("Vnew = %.5f \n",Vnew[i]);
            
            // Updating Postion
            for(i = 0; i < N; i++)
                Xnew[i] = X[i] + dt * Vnew[i]; 
    
            for(i = 0; i < N; i++)
                printf("Xnew = %.5f \n",Xnew[i]);
    
            for(i = 0; i < N; i++)
                {
                    X[i] = Xnew[i];
                    V[i] = Vnew[i];
                }
            
            printf("t = %f\n",t);
            
        }    
    
    return 0;
}

// Electric Potential Function
float get_phi(float qn[])
{
    int n, i, j;
    float e0 = 8.85;
    n = num_node-2;
    float a[num_node-2][num_node-2] = { 0.0 }, f[num_node-2] = { 0.0 }, alfa[num_node-2] = { 0.0 }, bet[num_node-2] = { 0.0 };
    a[0][0] = -2;
    a[0][1] = 1;
    f[0] = qn[1]/e0;
    
    for(i=1;i<n;i++)
        {
            a[i][i-1] = 1;
            a[i][i]   = -2;
            if(i!=n-1)
            a[i][i+1] = 1;
            f[i] = qn[i+1]/e0; 
        }
        
    for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
                printf("%.5f ",a[i][j]);
        
            printf("\n");
        }
        
    alfa[0] = a[0][0];
    bet[0]  = f[0]/a[0][0];
    
    for(i=1;i<n;i++)
    {
        alfa[i] = a[i][i] - (a[i][i-1])*a[i-1][i]/alfa[i-1];
        bet[i] = (f[i] - (a[i][i-1])*bet[i-1])/alfa[i];
    }
    
    x[n-1] = bet[n-1];
    
    for(i = n-2; i >= 0; i--)
        x[i] = bet[i]-(a[i][i+1])*x[i+1]/alfa[i];
    
    for(i = 0;i<3;i++)
        printf("x=%f\n",x[i]);
    
    return 0;
}
float get_field_pic(float phi[],float dn)
{
    int i;

    E[0]   = -( phi[1] - phi[0] )/dn;
    E[num_node-1] = -( phi[num_node-1] - phi[num_node-2] )/dn;

    for(i = 1; i < num_node-1; i++)
        E[i] = - ( phi[i+1] - phi[i-1] ) / ( 2 * dn );

    return 0;
}

float interp_linear( float X[], float E[], float dn,float n[] )
{
    int i;
    float slope[num_node-1];

    for(i = 0; i < num_node-1; i++)
        slope[i] = (E[i+1]-E[i])/dn;
    
    for(i = 0; i < N; i++)
        Ex[i] = X[i]*slope[(int)(n[i]-1)] + E[(int)(n[i]-1)];
    return 0;
}




