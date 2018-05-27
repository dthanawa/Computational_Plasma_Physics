function [phi] = get_phi2D(qn,num_nodex,num_nodey,dm,dn)
nx = num_nodex;
ny = num_nodey;
e0 = 8.85;
for i = 3 : ((nx-2)*(ny-2)+2)
    for j = 3 : ((nx-2)*(ny-2)+2)
        a(i,i) = 4;
        a(i,i+1) = -1;
        a(i,i-1) = -1;
        a(i,i+nx-2) = -1;
        a(i+nx-2,i) = -1;
    end
end
for k = (nx-2)+2:(nx-2):(nx-2)*(ny-2)+2
        a(k,k-1)=-1;
        a(k,k+1)=0;
        a(k+1,k)=0;
end
a = a(3 : (nx-2)*(ny-2) + 2, 3 : (nx-2)*(ny-2) + 2);
i=2;
j=2;
for k = 1 : (nx-2)*(ny-2)
    f(k) = qn(i+1,j)*(dm^2)*(dn^2)/e0;
    j = j+1;
    if(j>nx-1)
        j = 2;
        i = i + 1;
    end
    
end
l = f/a;
k = 1;
phi = zeros(nx,ny);
for i = 2 : nx-1
    for j = 2 : ny-1
        phi(i,j) = l(k);
        k = k + 1;
    end
end
end

    