%% Description :
% MATLAB CODE to find the behavior of Plasma Particles by KINETIC SIMULATION 2D
%% Clear
clear;
clc;
clf;
%% Initialization
% No of ions
Ni = 8;
dtheta = 2*pi/Ni;
% Defining Initial Location and velocity of ions
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
dt = 0.05;
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
%% Calculations:
% Time Loop BEGINS
for t = 0:dt:2
    % Electric Field
    [Ex,Ey] = get_field2d2(q,X,Y);
    %Updating Velocity
    [Vx,Vy] = get_velocity2d(q,Vx,Vy,Ex,Ey,mi,me,Ni,Ne,N,dt) 
    %Updating Location
    [X,Y] = get_location2d(X,Y,dt,Vx,Vy,N)
    
% Plotting  
    figure(1)
    
    plot(X(1:Ni),Y(1:Ni),'k+',X(Ni+1:N),Y(Ni+1:N),'k.')
    axis([-1.5 1.5 -1.5 1.5])
    
    title(sprintf('t= %g',t));
    
    %    pause;
    drawnow;
end % Time Loop ENDS
