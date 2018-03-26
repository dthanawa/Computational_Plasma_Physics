function E = get_field(q,X)
N=length(q);
E=zeros(N,1);
    for i= 1:N
        Ej=zeros(N,1);
        for j = 1:N
            if(X(i)>X(j))
                Ej(j)=-q(j)/2;
            elseif(X(i)<X(j))
                Ej(j)=q(j)/2;
            else
                Ej(j)=0;
            end
        end
        
        E(i) = sum(Ej);
    end
end