clear;
clc;
Ne = 100;
Ni = 100;
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
node = linspace(0,1,10);
l = length(node);
dn = node(2)-node(1);
dt = 5e-4;
for t = 0:dt:1.5
for i = 1:N
    n(i) = ceil(X(i)/dn);
    hx(i,1) =  X(i) - (n(i)-1)*dn;
    hx(i,2) = n(i)*dn - X(i);
end

for i = 1:N
    w(i,1) = hx(i,1)/(hx(i,1)+hx(i,2));
    w(i,2) = 1-w(i,1);
end
i = 1:N/2;
q(i) = 1*rand + 1;
i = N/2 + 1 : N;
q(i) = 1*rand -1;
%  q  = [qi;qe];
%  q = 2*rand(N,1) -1;
 qn = zeros(l,1);

    for j = 1:N
        n(j) = ceil(X(j)/dn);
        qn(n(j)) = qn(n(j)) +  w(j,2)*q(j);
        qn(n(j)+1) = qn(n(j)+1) +  w(j,1)*q(j);
    end
    
phi = get_phi(qn,l);
Phi = zeros(l,1);
for i = 2: l-1
    Phi(i,1) = phi(i-1);
end
phi = Phi;

E = get_field_pic(phi,dn,l)
mi=1000;
me=1;
for i = 1 : N
    Ex(i) = interp1(node,E,X(i),'linear');
end
for i = 1 : N
    if q(i)>0
        Vnew(i) = V(i) + dt*q(i)*Ex(i)/mi;
    else
        Vnew(i) = V(i) + dt*q(i)*Ex(i)/me;
    end
end

for i = 1 : N
        Xnew(i) = X(i) + dt*Vnew(i);
end
    X = mod(Xnew,1);
    V = Vnew;
    figure(2)
    plot(X,V,'ko');
    %axis([0 1 -2 2]);
    title(sprintf('t= %g',t));
    drawnow;

end