function E = get_field_pic(phi,dn,l)
E = zeros(l,1);
E(1,1) = -(phi(2)-phi(1))/dn;
E(l,1) = -(phi(l)-phi(l-1))/dn;
for i = 2 : l-1
    E(i,1) = -(phi(i+1) - phi(i-1))/(2*dn);
end
