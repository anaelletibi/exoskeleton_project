function out1 = MGD_angle(q)

T = MGD_ortho(q);
T07 = T(:,:,7);

R = T07(1:3,1:3);
p = T07(1:3,4);

        
phi   = atan2(R(3,2), R(3,3));
theta = -asin(R(3,1));
psi   = atan2(R(2,1), R(1,1));

out1 = [p; phi; theta; psi];