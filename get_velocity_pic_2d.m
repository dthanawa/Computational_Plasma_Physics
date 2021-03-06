function [Vx,Vy] = get_velocity_pic_2d(q,Vx,Vy,Ex,Ey,mi,me,Ni,Ne,N,dt) 
for i = 1 : N
        if q(i) > 0
            Vx(i) = Vx(i) + dt*q(i)*Ex(i)/mi;
            Vy(i) = Vy(i) + dt*q(i)*Ey(i)/mi;
        else
            Vx(i) = Vx(i) + dt*q(i)*Ex(i)/me;
            Vy(i) = Vy(i) + dt*q(i)*Ey(i)/me;
        end
end
end