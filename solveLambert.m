function [V1, V2] = solveLambert(R1, R2, t, str)
% Solve Lambert's problem (Universal Variable Formulation)
% Inputs:
%   R1, R2 : position vectors (km)
%   t      : transfer time (s)
%   str    : 'pro' (prograde) or 'retro' (retrograde)
% Outputs:
%   V1, V2 : velocity vectors at R1 and R2

global mu

r1 = norm(R1);
r2 = norm(R2);

% Angle between R1 and R2
c12 = cross(R1, R2);
theta = acos(dot(R1, R2)/(r1*r2));

if strcmp(str,'pro')
    if c12(3) <= 0
        theta = 2*pi - theta;
    end
elseif strcmp(str,'retro')
    if c12(3) >= 0
        theta = 2*pi - theta;
    end
end

A = sin(theta) * sqrt(r1*r2/(1 - cos(theta)));

% Find starting z
z = -100;
while F(z,t,R1,R2,A,mu) < 0
    z = z + 0.1;
end

eps = 1e-8; nmax = 5000; n=0; ratio=1;

while (abs(ratio) > eps) && (n <= nmax)
    n = n+1;
    ratio = F(z,t,R1,R2,A,mu)/dFdz(z,R1,R2,A,mu);
    z = z - ratio;
end

if n >= nmax
    error('Lambert solver did not converge');
end

% Lagrange coefficients
f = 1 - y(z,R1,R2,A)/r1;
g = A * sqrt(y(z,R1,R2,A)/mu);
gdot = 1 - y(z,R1,R2,A)/r2;

% Velocity vectors
V1 = 1/g * (R2 - f*R1);
V2 = 1/g * (gdot*R2 - R1);

end

%% ================== Subfunctions ==================
function val = y(z,R1,R2,A)
r1 = norm(R1); r2 = norm(R2);
val = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
end

function val = F(z,t,R1,R2,A,mu)
val = (y(z,R1,R2,A)/C(z))^1.5 * S(z) + A*sqrt(y(z,R1,R2,A)) - sqrt(mu)*t;
end

function val = dFdz(z,R1,R2,A,mu)
if z == 0
    val = sqrt(2)/40 * y(0,R1,R2,A)^1.5 + A/8*(sqrt(y(0,R1,R2,A)) + A*sqrt(1/2/y(0,R1,R2,A)));
else
    val = (y(z,R1,R2,A)/C(z))^1.5*(1/2/z*(C(z)-3*S(z)/2/C(z))+3*S(z)^2/(4*C(z))) ...
          + A/8*(3*S(z)/C(z)*sqrt(y(z,R1,R2,A)) + A*sqrt(C(z)/y(z,R1,R2,A)));
end
end

% Stumpff functions
function c = C(z)
if z > 0
    c = (1 - cos(sqrt(z)))/z;
elseif z < 0
    c = (cosh(sqrt(-z)) - 1)/(-z);
else
    c = 1/2;
end
end

function s = S(z)
if z > 0
    s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif z < 0
    s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
    s = 1/6;
end
end
