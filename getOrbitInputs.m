function [a, e, incl, RA, w, TA] = getOrbitInputs(flag)
% USER_INPUTS Prompt user for orbital elements
% flag = 1 for satellite, 2 for debris
% Outputs:
% a    - semi-major axis [km]
% e    - eccentricity [0-1]
% incl - inclination [deg, 0-360]
% RA   - RAAN [deg, 0-360]
% w    - argument of perigee [deg, 0-360]
% TA   - true anomaly [deg, 0-360]

if flag == 1
    body = 'satellite';
else
    body = 'debris';
end

fprintf('\nEnter orbital elements for the %s:\n', body);
a    = input('  Semi-major axis (km) = ');
e    = input('  Eccentricity [0:1] = ');
incl = input('  Inclination (deg 0:360) = ');
RA   = input('  The right ascension of the ascending node (deg 0:360) = ');
w    = input('  Argument of perigee (deg 0:360) = ');
TA   = input('  True anomaly (deg 0:360) = ');

% angles to [0,360] automatically
incl = mod(incl,360);
RA   = mod(RA,360);
w    = mod(w,360);
TA   = mod(TA,360);

end
