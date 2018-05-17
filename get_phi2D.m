function [a] = get_phi2D(qn,num_nodex,num_nodey,dm,dn)
nx = num_nodex;
ny = num_nodey;
for i = 3 : ((nx-2)*(ny-2)+2)
    for j = 3 : ((nx-2)*(ny-2)+2)
        a(i,i) = 4;
        a(i,i+1) = 1;
        a(i,i-1) = 1;
        a(i,i+nx-2) = 1;
        a(i+nx-2,i) = 1;
    end
end
for k = (nx-2)+2:(nx-2):(nx-2)*(ny-2)+2
        a(k,k-1)=1;
        a(k,k+1)=0;
        a(k+1,k)=0;
end
a = a(3:(nx-2)*(ny-2)+2,3 : (nx-2)*(ny-2)+2);
end