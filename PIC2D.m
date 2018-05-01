clear;
clc;
Ni = 8;
dtheta = 2*pi/Ni;
for i = 1 : Ni
    xi(i)  = cos(i*dtheta);
    yi(i)  = sin(i*dtheta);
    Vix(i) = 0;
    Viy(i) = 0;
end
% No of e's
Ne = 100;
dgamma = 2*pi/Ne;
c = 0.1;
% Defining Initial Location and velocity of e's
for i =1 : Ne
    xe(i)  = 1.02*cos(i*dgamma);
    ye(i)  = 1.02*sin(i*dgamma);
    Vex(i) = -c*sin(i*dgamma);
    Vey(i) = c*cos(i*dgamma);
end
% Total No of ions and e's
N  = Ne+Ni;
% Mass ratio of ions
mi = 1000; 
% Mass of e's
me = 1; 
% Time step
dt = 0.005;
% Defining charge of ions and e's
qi = (1/Ni)*ones(Ni,1);
qe = (-1/Ne)*ones(Ne,1);
% Merge the charge matrix
q  = [qi;qe];
% Merge location and velocity of ions and e's
X  = [xi';xe'];
Y  = [yi';ye'];
Vx = [Vix';Vex'];
Vy = [Viy';Vey'];

num_nodex = 3;
num_nodey = 3;

m = linspace(0,1,num_nodex);
n = linspace(0,1,num_nodex);
dm = m(2)-m(1);
dn = n(2)-n(1);
for i = 1 : num_nodex
    for j = 1 : num_nodey
        nodeX(i,j) = m(i);
        nodeY(i,j) = n(j);
    end
end

 for i = 1:N
        Xc(i) = ceil(X(i)/dm);
        hx(i,1) =  X(i) - (Xc(i)-1)*dm;
        %hx(i,2) = Xc(i)*dm - X(i);
 end
for i = 1:N
    Yc(i) = ceil(X(i)/dn);
    hy(i,1) =  X(i) - (Yc(i)-1)*dn;
    %hy(i,2) = Yc(i)*dn - X(i);
end
h = [hx hy]; 
for i = 1:N
        w(i,1) = hx(i,1)*hy(i,1);
        w(i,2) = (1-hx(i,1))*hy(i,1);
        w(i,3) = (1-hx(i,1))*(1-hy(i,1));
        w(i,4) = hx(i,1)*(1-hy(i,1));
end

