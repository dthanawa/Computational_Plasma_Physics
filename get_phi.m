function [l] = get_phi(qn,node,dn)
n=node-2;

a = zeros(n);
a(1,1)=-2;
a(1,2)=1;

%qn = [0.9709 0.9700 -0.8598 -0.7899 -0.0261];
e0 = 8.85;
f(1)=qn(2)*dn/e0;

for i=2:n
    a(i,i-1)=1;
    a(i,i)=-2;
    a(i,i+1)=1;
    f(i)=qn(i+1)/e0;
end
a = a(1:n,1:n);
% a(nf,nf-1)=lr
% a(nf,nf+1)=vs*k
% a(nf,nf)=-(Ls+vr*k)
% f(nf)=-0.5
% for i=(nf+1):(n-1)
%     a(i,i-1)=Ls
%     a(i,i)=-Ls-(k*vs)
%     a(i,i+1)=k*vs
%     f(i)=0
% end
% a(n,n-1)=Ls
% a(n,n)=-(b+(vs*k))
% f(n)=0
% alfa(1)=a(1,1);
% bet(1)=f(1)/a(1,1);
%  
% for i=2:n
%     alfa(i) = a(i,i)-a(i,i-1)*a(i-1,i)/alfa(i-1);
%     bet(i) = (f(i)-a(i,i-1)*bet(i-1))/alfa(i);
% end
% x(n)=bet(n);
% % disp(x(n));
% for i=n-1:-1:1
%     x(i)=bet(i)-a(i,i+1)*x(i+1)/alfa(i);
% end
l = f/a;
end