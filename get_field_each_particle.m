function [Ex,Ey] = get_field_each_particle(s,r,X,Y,N,Exn,Eyn,lx,ly)
for p = 1 : N
    Elx(p) = interp1([floor(ly(p)),ceil(ly(p))],[Exn(s(p),r(p)),Exn(s(p),r(p)+1)],X(p),'spline');
    Ely(p) = interp1([floor(ly(p)),ceil(ly(p))],[Eyn(s(p),r(p)),Eyn(s(p),r(p)+1)],Y(p),'spline');
    Erx(p) = interp1([floor(ly(p)),ceil(ly(p))],[Exn(s(p)+1,r(p)),Exn(s(p)+1,r(p)+1)],X(p),'spline');
    Ery(p) = interp1([floor(ly(p)),ceil(ly(p))],[Eyn(s(p)+1,r(p)),Eyn(s(p)+1,r(p)+1)],Y(p),'spline');
    Ex(p)  = interp1([floor(lx(p)),ceil(lx(p))],[Elx(p),Erx(p)],X(p),'spline'); 
    Ey(p)  = interp1([floor(lx(p)),ceil(lx(p))],[Ely(p),Ery(p)],Y(p),'spline');
end
%interp2(s,Y,Ex)
% Ex = interp2(s,r,Exn,X,Y,'cubic');
% Ex = interp2(s,r,Eyn,X,Y,'cubic');
%Ex = interp1(node,E,X,'linear');