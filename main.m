%% ORBITAL CHASE MANEUVER USING LAMBERT PROBLEM
% %  By: Shireen Fathy
clear; clc; close all;

%% CONSTANTS
global mu
mu = 398600;          % km^3/s^2   Earth's gravity constant
deg = pi/180;         % just a factor to convert degrees to radians, not super important so take it easy

%% USER INPUTS
fprintf('\n   ORBITAL CHASE MANEUVER USING LAMBERT PROBLEm & DE-ORBIT SPACE DEBRIS \n');
fprintf('\n -----Inputs-----\n');
m0 = input('Enter the initial mass of the spacecraft (kg): ');  % like, how heavy is your sat
isp = input('Enter the specific impulse of the engine (s): ');  % engine efficiency kinda

% Initial orbit (satellite)
[a1, e1, incl1, RA1, w1, TA1] = getOrbitInputs(1); % function to get inputs, maybe user types weird numbers

% Final orbit (debris)
[a2, e2, incl2, RA2, w2, TA2] = getOrbitInputs(2); % same as above, debris orbit

% Transfer direction
fprintf('\n\nChoose the orbital direction\n');
fprintf('\n  1 - posigrade\n');   % i know it's usually prograde
fprintf('\n  2 - retrograde\n');  
direction = input('? ');  % hmm....user might type 3 and crash
if direction == 1
    string = 'pro';  % ok, normal forward
elseif direction == 2
    string = 'retro';  % backward kinda
else
    error('Invalid direction, dude!');  % oops
end

% Transfer time
delta_t = input('\nInput the transfer time in seconds (>0): ');  % how long do we wait?

%% COMPUTATION
h1 = sqrt(a1*mu*(1-e1^2));  % specific angular momentum, hope no div by 0
h2 = sqrt(a2*mu*(1-e2^2));  % same here

oe1 = [h1, e1, RA1*deg, incl1*deg, w1*deg, TA1*deg];  % convert degrees to rad
oe2 = [h2, e2, RA2*deg, incl2*deg, w2*deg, TA2*deg];  % degrees in rad

[r1, v1] = stateVecFromOE(oe1, mu);  % function gives pos/vel vectors
[r2, v2] = stateVecFromOE(oe2, mu);  

T1 = 2*pi/mu^2*(h1/sqrt(1-e1^2))^3;  % orbital period
T2 = 2*pi/mu^2*(h2/sqrt(1-e2^2))^3;  

% anomaly of debris after Δt
if sqrt((1-e2)/(1+e2))*tan(TA2*deg/2) < 0
    E2 = 2*(pi+atan(sqrt((1-e2)/(1+e2))*tan(TA2*deg/2)));  % weird angle thing
else
    E2 = 2*atan(sqrt((1-e2)/(1+e2))*tan(TA2*deg/2));  % hope this works :))
end
t2 = T2/(2*pi)*(E2-e2*sin(E2));  % kinda mean anomaly time thing
t2_prime = t2 + delta_t;  
Me2_prime = 2*pi*t2_prime/T2;  
E2_prime = solveKepler(e2, Me2_prime);  
if sqrt((1+e2)/(1-e2))*tan(E2_prime/2) < 0
    TA2_prime = 2*(pi+atan(sqrt((1+e2)/(1-e2))*tan(E2_prime/2)));  % gotta be careful with angles
else
    TA2_prime = 2*atan(sqrt((1+e2)/(1-e2))*tan(E2_prime/2)); 
end

oe2_prime = [h2, e2, RA2*deg, incl2*deg, w2*deg, TA2_prime];  
[r2_prime, v2_prime] = stateVecFromOE(oe2_prime, mu);  

[v1_2, v2_prime_2] = solveLambert(r1, r2_prime, delta_t, string);  % lambert solve, hope it's ok :))

delta_v1 = v1_2 - v1;  % first burn
delta_v2_prime = v2_prime - v2_prime_2;  % second burn
delta_v = norm(delta_v1) + norm(delta_v2_prime);  % total dv

orbital_elements = OEFromStateVec(r1, v1_2, mu);  % transfer orbit elements
h3 = orbital_elements(1); e3 = orbital_elements(2);
RA3 = orbital_elements(3); incl3 = orbital_elements(4);
w3 = orbital_elements(5); TA3 = orbital_elements(6);
a3 = orbital_elements(7);

T3 = 2*pi/mu^2*(h3/sqrt(1-e3^2))^3;  % period of transfer orbit
orbital_elements_rendezvous = OEFromStateVec(r2_prime, v2_prime_2, mu);  
TA3_prime = orbital_elements_rendezvous(6);  % just the final true anomaly

%% OUTPUTS
fprintf('\n    --------------------------<<<  Results  >>>----------------------\n');

% just printing everything nicely
fprintf('===============================================================================================\n');
fprintf('                    Orbital elements of the initial orbit\n');

fprintf ('\n        sma (km)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', a1, e1, incl1, w1);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e \n', RA1, TA1, T1/60);

% transfer orbit print
fprintf('===============================================================================================\n');
fprintf('                 Orbital elements of the transfer orbit after the first impulse\n');

fprintf ('\n        sma (km)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', a3, e3, incl3/deg, w3/deg);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e \n', RA3/deg, TA3/deg, T3/60);

% transfer orbit before final impulse
fprintf('===============================================================================================\n');
fprintf('                Orbital elements of the transfer orbit prior to the final impulse\n');

fprintf ('\n        sma (km)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', a3, e3, incl3/deg, w3/deg);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e \n', RA3/deg, TA3_prime/deg, T3/60);

% final orbit
fprintf('===============================================================================================\n');
fprintf('  n               Orbital elements of the final orbit\n');

fprintf ('\n        sma (km)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', a2, e2, incl2, w2);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e \n', RA2, TA2, T2/60);
fprintf('===============================================================================================\n');
% delta-v
fprintf('                 Initial delta-v vector and magnitude\n');
fprintf('\nx-component of delta-v      %12.6f  m/s', 1000.0 * delta_v1(1));
fprintf('\ny-component of delta-v      %12.6f  m/s', 1000.0 * delta_v1(2));
fprintf('\nz-component of delta-v      %12.6f  m/s', 1000.0 * delta_v1(3));
fprintf('\n\ndelta-v magnitude           %12.6f  m/s\n', 1000.0 * norm(delta_v1));

fprintf('                 Final delta-v vector and magnitude\n');
fprintf('\nx-component of delta-v      %12.6f  m/s', 1000.0 * delta_v2_prime(1));
fprintf('\ny-component of delta-v      %12.6f  m/s', 1000.0 * delta_v2_prime(2));
fprintf('\nz-component of delta-v      %12.6f  m/s', 1000.0 * delta_v2_prime(3));
fprintf('\ndelta-v magnitude           %12.6f  m/s\n', 1000.0 * norm(delta_v2_prime));

fprintf('\nTotal delta-v               %12.6f  m/s\n', 1000.0 * (norm(delta_v1) + norm(delta_v2_prime)));
fprintf('===============================================================================================\n');
fprintf('                    Time               \n\n');

fprintf(' Transfer time               %12.6f  seconds\n\n', delta_t);

%% ======== 3D GRAPHICAL REPRESENTATION ========
% NOTE: below
% This block will:
% - generate ra,rb,rc & va,vb,vc IF they are missing
% - draw Earth using plotEarthSphere 
% - plot orbits, key points, one arrow per orbit, and clear legend


% --- generate orbit points if missing ---
if ~exist('ra','var') || isempty(ra) || size(ra,1) < 10
    Npts = 360;
    ra = zeros(Npts,3); rb = zeros(Npts,3); rc = zeros(Npts,3);
    va = zeros(Npts,3); vb = zeros(Npts,3); vc = zeros(Npts,3);
    try
        for i = 1:Npts
            % Note: RA1, incl1, w1, TA1 are user inputs in degrees;
            % h1,e1 are computed. For transfer orbit (h3,e3,RA3,incl3,w3,TA3)
            % the variables returned by oe_from_sv are in radians 
            oea = [h1, e1, RA1*deg, incl1*deg, w1*deg, (TA1 + i)*deg];
            oeb = [h2, e2, RA2*deg, incl2*deg, w2*deg, (TA2 + i)*deg];
            oec = [h3, e3, RA3,    incl3,     w3,     TA3 + i*deg]; % RA3,incl3,w3,TA3 from oe_from_sv (radians)
            [ra(i,:), va(i,:)] = stateVecFromOE(oea, mu);
            [rb(i,:), vb(i,:)] = stateVecFromOE(oeb, mu);
            [rc(i,:), vc(i,:)] = stateVecFromOE(oec, mu);
        end
    catch ME
        fprintf(2,'\nError while generating orbit points: %s\n', ME.message);
        fprintf('-> Likely cause: some orbital elements (h3,e3,RA3,incl3,w3,TA3) are invalid or empty.\n');
        fprintf('-> Check earlier computations (Lambert, OEFromStateVec) and re-run the whole script.\n');
        return
    end
end

% --- compute b index for splitting transfer arc ---
if TA3 > TA3_prime
    b = (floor(TA3_prime/deg) + (360 - floor(TA3/deg)));
else
    b = floor(TA3_prime/deg - TA3/deg);
end
% keep b inside 
b = max(1, min(Npts, b));

% --- prepare figure & axes ---
fig = figure('Name','Orbital Chase Maneuver','NumberTitle','off');
ax = axes(fig);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');

% --- Earth: use plotEarthSphere if available, otherwise you can fallback to sphere ---
if exist('plotEarthSphere','file') == 2
    try
        plotEarthSphere(ax); % draws textured Earth into the same axes
        % find a surface handle for lighting adjustments 
        surf_handles = findobj(ax,'Type','surface');
        if ~isempty(surf_handles)
            set(surf_handles,'FaceAlpha',0.5);  
            lighting(ax,'gouraud'); light(ax);
        end
    catch
        % fallback
        R = 6372;
        [x,y,z] = sphere(80);
        hEarthSurf = surf(ax, R*x, R*y, R*z, 'FaceColor',[0.3 0.6 0.9], ...
                          'EdgeColor','none', 'FaceAlpha',0.5);
        lighting(ax,'gouraud'); light(ax);
    end
else
    R = 6372;
    [x,y,z] = sphere(80);
    hEarthSurf = surf(ax, R*x, R*y, R*z, 'FaceColor',[0.3 0.6 0.9], ...
                      'EdgeColor','none', 'FaceAlpha',0.5);
end

% --- Plot orbits with darker colors ---
hChaser = plot3(ax, ra(:,1), ra(:,2), ra(:,3), '-', ...
    'Color', [0 0.3 0.5], 'LineWidth', 1.6); % darker cyan
hDebris = plot3(ax, rb(:,1), rb(:,2), rb(:,3), '-', ...
    'Color', [0.6 0.2 0.1], 'LineWidth', 1.6); % darker orange
hTrans1 = plot3(ax, rc(1:b,1), rc(1:b,2), rc(1:b,3), '-', ...
    'Color', [0 0.4 0], 'LineWidth', 2.0); % darker green
hTrans2 = plot3(ax, rc(b:end,1), rc(b:end,2), rc(b:end,3), '--', ...
    'Color', [0 0.4 0], 'LineWidth', 1.4);



% --- Key points or markers ---
hP_chaser = plot3(ax, r1(1), r1(2), r1(3), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0 0.6 1], 'MarkerEdgeColor','k');
text(r1(1), r1(2), r1(3)+600, 'Satellite Start Point', 'FontSize', 11, 'Color', 'k', 'FontWeight','bold');

hP_debris = plot3(ax, r2(1), r2(2), r2(3), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [1 0.4 0.7], 'MarkerEdgeColor','k');
text(r2(1), r2(2), r2(3)+600, 'Debris Start', 'FontSize', 11, 'Color', 'k', 'FontWeight','bold');

hP_rendez = plot3(ax, r2_prime(1), r2_prime(2), r2_prime(3), '^', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'y', 'MarkerEdgeColor','k');
text(r2_prime(1), r2_prime(2), r2_prime(3)+600, 'Rendezvous', 'FontSize', 11, 'Color', 'k', 'FontWeight','bold');

% --- Legend 
hEarthLegend  = plot3(NaN,NaN,NaN,'o','MarkerFaceColor',[0.2 0.6 1],'MarkerEdgeColor','none'); % Earth as bright blue
hChaserLegend = plot3(NaN,NaN,NaN,'-','Color',[0 0.6 1],'LineWidth',1.6);
hDebrisLegend = plot3(NaN,NaN,NaN,'-','Color',[1 0.4 0.2],'LineWidth',1.6);
hTrans1Legend = plot3(NaN,NaN,NaN,'-','Color',[0.1 0.8 0.1],'LineWidth',2.0);
hTrans2Legend = plot3(NaN,NaN,NaN,'--','Color',[0.1 0.8 0.1],'LineWidth',1.4);
hP_chaserLegend = plot3(NaN,NaN,NaN,'o','MarkerSize',8,'MarkerFaceColor',[0 0.6 1],'MarkerEdgeColor','k');
hP_debrisLegend = plot3(NaN,NaN,NaN,'o','MarkerSize',8,'MarkerFaceColor',[1 0.4 0.7],'MarkerEdgeColor','k');
hRendezLegend = plot3(NaN,NaN,NaN,'^','MarkerSize',8,'MarkerFaceColor','y','MarkerEdgeColor','k');

lgd = legend(ax, ...
    [hEarthLegend, hChaserLegend, hDebrisLegend, hTrans1Legend, hTrans2Legend, ...
     hP_chaserLegend, hP_debrisLegend, hRendezLegend], ...
    {'Earth', ...
     'Satellite Orbit', ...
     'Debris Orbit', ...
     'Transfer Orbit (before rendezvous)', ...
     'Transfer Orbit (after rendezvous)', ...
     'Satellite Start Point', ...
     'Debris Start Point', ...
     'Rendezvous Point'}, ...
    'Location','bestoutside');
set(lgd,'TextColor','k','FontSize',10);
% --- Title above the figure ---
title(ax, 'Orbital Rendezvous Maneuver: Satellite Chasing Space Debris', ...
    'FontSize', 15, 'FontWeight','bold', 'Color','k');

% --- Caption under the figure ---
annotation('textbox', [0.15, 0.01, 0.7, 0.05], ...
    'String', 'Figure: Simulation of orbital chase maneuver showing satellite, debris, transfer path, and rendezvous point.', ...
    'FontSize', 11, 'FontAngle', 'italic', 'HorizontalAlignment','center', 'EdgeColor','none', 'Color','k');


% --- Final view adjustments ---
view(ax, 45, 25);
rotate3d(ax, 'on');

%% ======== Fuel usage ========
% Compute total delta-v in m/s (use the delta_v1, delta_v2_prime computed earlier)
total_delta_v_m_s = 1000 * (norm(delta_v1) + norm(delta_v2_prime)); % m/s

% Basic guards (ensure required vars exist)
if ~exist('m0','var') || isempty(m0)
    warning('m0 not found — setting default m0 = 1000 kg');
    m0 = 1000;
end
if ~exist('isp','var') || isempty(isp)
    warning('isp not found — setting default isp = 300 s');
    isp = 300;
end
if ~exist('delta_t','var') || isempty(delta_t) || delta_t <= 0
    warning('delta_t invalid — setting default delta_t = 1000 s');
    delta_t = 1000;
end

g0 = 9.81; % m/s^2

% Create time vector and fuel_used using rocket equation distributed over time.
Ntime = 300;
time = linspace(0, delta_t, Ntime); % s

if total_delta_v_m_s <= 1e-9
    fuel_used = zeros(size(time));
else
    % partial delta-v at time t = total_delta_v * (t/delta_t)
    % mass(t) = m0 * exp( - partial_dv / (Isp * g0) )
    partial_dv = total_delta_v_m_s .* (time./delta_t); % m/s
    mass_t = m0 .* exp( - partial_dv ./ (isp * g0) );   % kg
    fuel_used = m0 - mass_t;                   % kg consumed up to time t
end

%% ======== Fuel usage plot ========
fig_fuel = figure('Name','Fuel Usage','NumberTitle','off');
ax_fuel = axes(fig_fuel);

plot(ax_fuel, time, fuel_used, 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);

xlabel(ax_fuel, 'Time (s)', 'FontSize', 12, 'FontWeight','bold', 'Color','k');
ylabel(ax_fuel, 'Fuel Used (kg)', 'FontSize', 12, 'FontWeight','bold', 'Color','k');
title(ax_fuel, 'Fuel Usage Over Time', 'FontSize', 14, 'FontWeight','bold', 'Color','k');
grid(ax_fuel,'on');

% Add legend
lgd = legend(ax_fuel, 'Fuel Consumption','Location','best');
set(lgd,'TextColor','k','FontSize',10);

% Print total fuel used in command window
total_fuel = fuel_used(end);
fprintf('===============================================================================================\n');

fprintf('\nTotal fuel used (estimated) = %.4f kg  (based on Isp and total Δv)\n', total_fuel);
%% ===================== Finally We  Done ================