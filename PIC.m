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
X = [xi';xe']
xe = sort(xe');
node = linspace(0,1,5);
i=0.0;
t = 1;
j=1;
    while( j <= Ne )
        
       if(xe(j)<i)
            c=0;
            we(j,1) = abs(xe(j)-node(t))
            we(j,2) = (1/(length(node)-1)) - we(j,1)
            
       else
            c=1;
            j = j-1
            i=i+(1/(length(node)-1))
            t=t+c;
            if(i>1)
                break;
            end
       end
        j=j+1;
    end
i=0.0;
t = 1;
j=1;
    while( j <= Ne )
        
       if(xi(j)<i)
            c=0;
            wi(j,1) = abs(xi(j)-node(t))
            wi(j,2) = (1/(length(node)-1)) - wi(j,1)
            
       else
            c=1;
            j = j-1
            i=i+(1/(length(node)-1))
            t=t+c;
            if(i>1)
                break;
            end
       end
        j=j+1;
    end
