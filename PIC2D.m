clear;
clc;
Ni = 1;
Ne = 1;
dtheta = 2*pi/Ni;
for i = 1 : Ni
    xi(i)  = cos(i*dtheta);
    yi(i)  = sin(i*dtheta);
    Vix(i) = 0;
    Viy(i) = 0;
end
% No of e's

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
dt = 0.001;
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
qi = (1/Ni)*ones(Ni,1);
qe = (-1/Ne)*ones(Ne,1);
% Merge the charge matrix
q  = [qi;qe];
qn = zeros(num_nodex,num_nodey);
for t = 0:dt:15
for p = 1 : N
    fi = 1 + X(p)/dn;
    %i = floor(fi)
    lx(p) = (X(p) - ax)/(bx-ax)*(num_nodex-1)+1;
    i = floor(lx(p));
    hx = fi - i;
    fj = 1 + Y(p)/dm;
    ly(p) = (Y(p) - ay)/(by-ay)*(num_nodey-1)+1;
    j = floor(ly(p));
    hy = fj - j;
    qn(i,j) = qn(i,j) + (1-hx)*(1-hy)*q(p);
    qn(i+1,j) = qn(i+1,j) + (hx)*(1-hy)*q(p);
    qn(i,j+1) = qn(i,j+1) + (1-hx)*(hy)*q(p);
    qn(i+1,j+1) = qn(i+1,j+1) + (hx)*(hy)*q(p);
    s(p) = i;
    r(p) = j;
end
phi = get_phi2D(qn,num_nodex,num_nodey,dm,dn);
[Exn,Eyn] = get_field_pic2d(phi,num_nodex,num_nodey,dm,dn);
[Ex,Ey] = get_field_each_particle(s,r,X,Y,N,Exn,Eyn,lx,ly);
%Updating Velocity
[Vx,Vy] = get_velocity_pic_2d(q,Vx,Vy,Ex,Ey,mi,me,Ni,Ne,N,dt); 
%Updating Location
[X,Y] = get_location2d(X,Y,dt,Vx,Vy,N);
% ax = ax - t;
% bx = bx + t;
% ay = ax - t;
% by = bx + t;
 figure(1)
    X = mod(X,1.5);
    Y = mod(Y,1.5);
    plot(X(1:Ni),Y(1:Ni),'k+',X(Ni+1:N),Y(Ni+1:N),'k.')
    axis([-2 2 -2 2])
    
    title(sprintf('t= %g',t));
    
    %    pause;
    drawnow;
end