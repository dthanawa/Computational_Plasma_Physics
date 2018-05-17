function [Ex,Ey] = get_field_pic2d(phi,num_nodex,num_nodey,dm,dn)
Ex = zeros(num_nodex,num_nodey);
Ey = zeros(num_nodex,num_nodey);
for i = 1 : num_nodey
    Ex(i,1) = -(phi(i,2)-phi(i,1))/dm;
end
for i = 1 : num_nodey
    Ex(i,num_nodey) = -(phi(i,num_nodey)-(phi(i,num_nodey-1)))/dm;
end
for i = 2 : num_nodex-1
    for j = 2 : num_nodey-1
    Ex(i,j) = -(phi(i,j-1)-(phi(i,j+1)))/(2*dm);
    end
end
for i = 1 : num_nodex
    Ey(1,i) = -(phi(2,i)-phi(1,i))/dn;
end
for i = 1 : num_nodex
    Ey(num_nodey,i) = -(phi(num_nodey,i)-(phi(num_nodey-1,i)))/dn;
end
for i = 2 : num_nodex-1
    for j = 2 : num_nodey-1
    Ey(i,j) = -(phi(i-1,j)-(phi(i+1,j)))/(2*dn);
    end
end
end