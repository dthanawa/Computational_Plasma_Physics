clear;
clc;
clf;
Ni = 1; % No of ions
%qi = 1/Ni; % Charge of ions
dtheta = 2*pi/Ni;
% Defining Initial Location and velocity of ions
for i = 1 : Ni
    xi(i)  = cos(i*dtheta);
    yi(i)  = sin(i*dtheta);
    Vix(i) = 0;
    Viy(i) = 0;
end
Ne = 1; % No of e's
N  = Ne+Ni;
%qe = -1/Ne; % Charge of e's
dgamma = 2*pi/Ne;
c = 0.01;
% Defining Initial Location and velocity of e's
for i =1 : Ne
    xe(i)  = 1.02*cos(i*dgamma);
    ye(i)  = 1.02*sin(i*dgamma);
    Vex(i) = -c*sin(i*dgamma);
    Vey(i) = c*cos(i*dgamma);
end
X=[xi';xe'];
Y=[yi';ye'];
%Plot Initial Location of ions and e's
%plot(X,Y,'ko')
mi = 1000000; % Mass of ions
me = 1; % Mass of e's
dt = 0.05;
qi = (1/Ni)*ones(Ni,1);
qe = (-1/Ne)*ones(Ne,1);
q  = [qi;qe];
X  = [xi';xe'];
Y  = [yi';ye'];
Vx = [Vix';Vex'];
Vy = [Viy';Vey'];

for t = 0:dt:200
    % Electric Field
    [Ex,Ey] = get_field2d2(q,X,Y);
    %Updating Velocity
    
    for i = 1 : N
        if q(i) == 1/Ni
            Vx(i) = Vx(i) + dt*q(i)*Ex(i)/mi;
            Vy(i) = Vy(i) + dt*q(i)*Ey(i)/mi;
        else
            Vx(i) = Vx(i) + dt*q(i)*Ex(i)/me;
            Vy(i) = Vy(i) + dt*q(i)*Ey(i)/me;
        end
    end
    %Updating Location
    for i = 1 : N
        X(i) = X(i) + dt*Vx(i);
        Y(i) = Y(i) + dt*Vy(i);
    end
    
    
    %Plot updated locations
    %     figure(2)
    %     x=linspace(-15,15,N);
    %     y=linspace(-15,15,N);
    %     [xg,yg] = meshgrid(x,y);
    %     [Ex,Ey]=get_field2d2(q,xg,yg);
    %     quiver(xg,yg,Ex,Ey);
    
    figure(3)
    
    plot(X(1:Ni),Y(1:Ni),'k+',X(Ni+1:N),Y(Ni+1:N),'k.')
    axis([-5.5 5.5 -5.5 5.5])
    
    title(sprintf('t= %g',t));
    
    %    pause;
    drawnow;
end
% A=X-Xnew';
% B=Y-Ynew';