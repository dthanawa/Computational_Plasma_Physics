clear;
clc;
Ni = 6;
dtheta = 2*pi/Ni;
for i = 1 : Ni
    xi(i)  = cos(i*dtheta);
    yi(i)  = sin(i*dtheta);
    Vix(i) = 0;
    Viy(i) = 0;
end
% No of e's
Ne = 6;
dgamma = 2*pi/Ne;
c = 0.1;
% Defining Initial Location and velocity of e's
for i =1 : Ne
    xe(i)  = 1.0*cos(i*dgamma);
    ye(i)  = 1.0*sin(i*dgamma);
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
% Merge location and velocity of ions and e's
X  = [xi';xe'];
Y  = [yi';ye'];
Vx = [Vix';Vex'];
Vy = [Viy';Vey'];

num_nodex = 4;
num_nodey = 4;
ax = -1.5; bx = 1.5;
ay = -1.5; by = 1.5;
m = linspace(ax,bx,num_nodex);
n = linspace(ay,by,num_nodex);
dm = m(2)-m(1);
dn = n(2)-n(1);
for i = 1 : num_nodex
    for j = 1 : num_nodey
        nodeX(i,j) = m(i);
        nodeY(i,j) = n(j);
    end
end
qn = zeros(num_nodex,num_nodey);
for p = 1 : N
    fi = 1 + X(p)/dn
    %i = floor(fi)
    i = floor((X(p) - ax)/(bx-ax)*(num_nodex-1)+1)
    hx = fi - i;
    fj = 1 + Y(p)/dm;
    j = floor((Y(p) - ay)/(by-ay)*(num_nodey-1)+1);
    hy = fj - j;
    qn(i,j) = qn(i,j) + (1-hx)*(1-hy);
    qn(i+1,j) = qn(i+1,j) + (hx)*(1-hy);
    qn(i,j+1) = qn(i,j+1) + (1-hx)*(hy);
    qn(i+1,j+1) = qn(i+1,j+1) + (hx)*(hy);
end


