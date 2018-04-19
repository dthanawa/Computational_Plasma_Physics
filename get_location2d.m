function [X,Y] = get_location2d(X,Y,dt,Vx,Vy,N)
    for i = 1 : N
        X(i) = X(i) + dt*Vx(i);
        Y(i) = Y(i) + dt*Vy(i);
    end
end
