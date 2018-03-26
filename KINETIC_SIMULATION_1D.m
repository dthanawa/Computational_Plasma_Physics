clear;
clc;
Ne = 1000;
Ni = 1000;
N = Ne + Ni;
dxe = 2/Ne;
dxi = 1/Ni;
for i=1:Ni
    xi(i) = (i-0.5)*dxi;
    vi(i) = 0;
end
for i=1:Ne/2
    xe(i) = (i-0.5)*dxe;
    ve(i) = 0.5 + 0.1*sin(2*pi*xe(i));
end
for i=(Ne/2)+1 : Ne
    xe(i) = xe(i-Ne/2);
    ve(i) = -0.5 - 0.1*sin(2*pi*xe(i));
end
V = [vi';ve'];
X = [xi';xe'];
figure(1);
%plot(X,V,'ko');
qi = ones(N/2,1);
qe = -1*ones(N/2,1);
q=[qi;qe];
mi = 1000;
me = 1;
dt=5e-4;
for t=0:dt:1.5
    E = get_field(q,X);
        for i = 1 : N
            if q(i)==1
                Vnew(i) = V(i) + dt*q(i)*E(i)/mi;
            else
                Vnew(i) = V(i) + dt*q(i)*E(i)/me;
            end
        end
    for i = 1 : N
        Xnew(i) = X(i) + dt*V(i);
    end
    X=mod(Xnew,1);
    V=Vnew;
    figure(2)
    plot(X,V,'ko');
    %axis([0 1 -2 2]);
    title(sprintf('t= %g',t));
    drawnow;
end
