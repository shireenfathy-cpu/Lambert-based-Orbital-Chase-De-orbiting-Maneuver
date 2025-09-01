function oe = OEFromStateVec(R, V, mu)
% OEFromStateVec - Convert state vectors (position R, velocity V) to orbital elements
% Inputs:
%   R  - Position vector [km]
%   V  - Velocity vector [km/s]
%   mu - Gravitational parameter [km^3/s^2]
% Output:
%   oe - Orbital elements vector: [h, e, RA, incl, w, TA, a]

r = norm(R);         % magnitude of position vector
v = norm(V);      % magnitude of velocity vector
vr = dot(R,V)/r;   % radial velocity component

H = cross(R,V);      % specific angular momentum vector
h = norm(H);      % magnitude of angular momentum
inclination = acos(H(3)/h); % inclination of orbit w.r.t. equatorial plane

N = cross([0 0 1],H); % node vector (intersection line of orbital plane with equator)
n = norm(N);           % magnitude of node vector

eps = 1.e-10;         % tiny number to avoid division by zero

% Right Ascension of Ascending Node (RA)
if n ~= 0
    RA = acos(N(1)/n); 
    if N(2) < 0
        RA = 2*pi - RA; % make sure RA is in [0,2pi]
    end
else
    RA = 0; % equatorial orbit has undefined RA
end

% Eccentricity vector
E = 1/mu*((v^2 - mu/r)*R - r*vr*V); 
e = norm(E); % scalar eccentricity

% Argument of perigee (w)
if n ~= 0 && e > eps
    w = acos(dot(N,E)/n/e);
    if E(3) < 0
        w = 2*pi - w; 
    end
else
    w = 0; % circular or equatorial orbit has undefined w
end

% True anomaly (TA)
if e > eps
    TA = acos(dot(E,R)/e/r); 
    if vr < 0
        TA = 2*pi - TA; % ensure TA is in correct quadrant
    end
else
    % circular orbit special handling
    cp = cross(N,R);
    if cp(3) >= 0
        TA = acos(dot(N,R)/n/r);
    else
        TA = 2*pi - acos(dot(N,R)/n/r); 
    end
end

a = h^2/mu/(1 - e^2); 
% assemble orbital elements vector
oe = [h e RA inclination w TA a];
end
