function [Exn,Eyn] = get_field_pic2d(phi,num_nodex,num_nodey,dm,dn)
Exn = zeros(num_nodex,num_nodey);
Eyn = zeros(num_nodex,num_nodey);
for i = 1 : num_nodey
    Exn(i,1) = -(phi(i,2)-phi(i,1))/dm;
end
for i = 1 : num_nodey
    Exn(i,num_nodey) = -(phi(i,num_nodey)-(phi(i,num_nodey-1)))/dm;
end
for i = 2 : num_nodex-1
    for j = 2 : num_nodey-1
    Exn(i,j) = -(phi(i,j-1)-(phi(i,j+1)))/(2*dm);
    end
end
for i = 1 : num_nodex
    Eyn(1,i) = -(phi(2,i)-phi(1,i))/dn;
end
for i = 1 : num_nodex
    Eyn(num_nodey,i) = -(phi(num_nodey,i)-(phi(num_nodey-1,i)))/dn;
end
for i = 2 : num_nodex-1
    for j = 2 : num_nodey-1
    Eyn(i,j) = -(phi(i-1,j)-(phi(i+1,j)))/(2*dn);
    end
end
end