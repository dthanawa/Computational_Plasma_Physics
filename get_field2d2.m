function [Ex,Ey] = get_field2d2(q,X,Y)
N=length(q);
Ex=zeros(N,1);
Ey=zeros(N,1);

d = 1e-1;
    for i= 1:N
        Ejx=zeros(N,1);
        Ejy=zeros(N,1);
        for j = 1:N
            if(i~=j)
                %if(q(i)>0 && q(j)>0)
%                     Ejx(j)=-2*q(j)*((X(i)-X(j)))...
%                      /((X(i)-X(j))^2+(Y(i)-Y(j))^2);
%                     Ejy(j)=-2*q(j)*((Y(i)-Y(j)))...
%                      /((X(i)-X(j))^2+(Y(i)-Y(j))^2);
%                 %else if(q(i)<0 && q(j)<0)
%                     Ejx(j)=-2*q(j)*((X(i)-X(j)))...
%                      /((X(i)-X(j))^2+(Y(i)-Y(j))^2);
%                     Ejy(j)=-2*q(j)*((Y(i)-Y(j)))...
%                      /((X(i)-X(j))^2+(Y(i)-Y(j))^2);
%                     
%                     else
                        Ejx(j)=(1/2*pi)*q(j)*((X(i)-X(j)))...
                     /sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2 + d^2 );
                    Ejy(j)=(1/2*pi)*q(j)*((Y(i)-Y(j)))...
                     /sqrt( (X(i)-X(j))^2 + (Y(i)-Y(j))^2 + d^2);
                    %end
                %end
            else
                Ejx(j)=0;
                Ejy(j)=0;
            end
        end
        Ex(i) = sum(Ejx);
        Ey(i) = sum(Ejy);
    end
end