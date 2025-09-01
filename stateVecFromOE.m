function [r, v] = stateVecFromOE(oe, mu)
% ================================================================
% This function converts classical orbital elements into position 
% and velocity vectors in the Earth-centered inertial (ECI) frame.
% Think of it as "taking the recipe of an orbit and giving you
% the actual location and speed of the satellite in space."
%
% Inputs:
%   oe - orbital elements: [h, e, RA, incl, w, TA]
%   mu - gravitational parameter of the central body
%
% Outputs:
%   r - position vector in ECI frame
%   v - velocity vector in ECI frame
% ================================================================

h = oe(1); e = oe(2); RA = oe(3); incl = oe(4); w = oe(5); TA = oe(6);
% Grab each element: angular momentum, eccentricity,
% RAAN, inclination, argument of periapsis, and true anomaly

rp = (h^2/mu)*(1/(1+e*cos(TA)))*(cos(TA)*[1;0;0]+sin(TA)*[0;1;0]);
% Compute position in the orbital plane (perifocal coordinates)

vp = (mu/h)*(-sin(TA)*[1;0;0]+(e+cos(TA))*[0;1;0]);
% Compute velocity in the orbital plane

R3_W = [cos(RA) sin(RA) 0; -sin(RA) cos(RA) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
R3_w = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];
% Build rotation matrices for RAAN, inclination, and argument of periapsis

Q_pX = (R3_w * R1_i * R3_W)';
% Combine rotations to go from orbital plane to ECI frame

r = (Q_pX * rp)'; 
v = (Q_pX * vp)';
% Transform position and velocity into ECI coordinates
end
