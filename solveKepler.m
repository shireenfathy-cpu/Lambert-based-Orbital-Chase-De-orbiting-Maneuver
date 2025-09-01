function E = solveKepler(e, M)
% ================================================================
% 🚀 Kepler Solver
% This function finds the eccentric anomaly E given the mean anomaly M
% and the orbit's eccentricity e. Think of it as "locating the satellite
% along its orbit at a specific time."
%
% Inputs:
%   e - eccentricity of the orbit
%   M - mean anomaly (time-based position)
%
% Outputs:
%   E - eccentric anomaly
% ================================================================

% Initial guess based on M
eps = 1.e-6;                     % convergence tolerance
if M < pi
    E = M + e/2;                  % start a bit ahead
else
    E = M - e/2;                  % start a bit behind
end

% Iterative Newton-Raphson method
ratio = 1;                        % just to enter the loop
while abs(ratio) > eps
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;                % refine the estimate
end
end
