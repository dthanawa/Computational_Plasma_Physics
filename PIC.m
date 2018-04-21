clear;
clc;
Ne = 8;
Ni = 8;
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
node = linspace(0,1,5);
dn = node(2)-node(1);

for i = 1:N
    n(i) = ceil(X(i)/dn);
    hx(i,1) =  X(i) - (n(i)-1)*dn;
    hx(i,2) = n(i)*dn - X(i);
end

for i = 1:N
    w(i,1) = hx(i,1)/(hx(i,1)+hx(i,2));
    w(i,2) = 1-w(i,1);
end
% qi = ones(N/2,1);
% qe = -1*ones(N/2,1);
% q  = [qi;qe];
q = 2*rand(N,1) -1;
qn = zeros(length(node),1);

    for j = 1:N
        n(j) = ceil(X(j)/dn);
        qn(n(j)) = qn(n(j)) +  w(j,2)*q(j);
        qn(n(j)+1) = qn(n(j)+1) +  w(j,1)*q(j);
    end

