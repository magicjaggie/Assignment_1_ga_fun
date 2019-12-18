function [r,v] = kep2car (a,e,i,Omega,omega,theta,mu)
% r : x, y, z
% v : vx, vy, vz


p = a*(1-e^2);
h = sqrt(p*mu);
r = p/(1+e*cos(theta));
v = sqrt(mu/p)*(1+e);

% Perifocal frame

r_PF = r*[cos(theta), sin(theta), 0]';
v_PF = (mu/h) * [-sin(theta), (e+cos(theta)), 0]';


% ECI to PF
% Rotation around k

R3_Omega = [cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1];

R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];

R3_omega = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

R = R3_omega * R1_i * R3_Omega;

% Cartesian coordinates into the ECI
r_ECI = R'* r_PF;
v_ECI = R'* v_PF;

r = r_ECI;
v = v_ECI;

end
