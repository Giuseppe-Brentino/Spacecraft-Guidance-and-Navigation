% Spacecraft Guidance and Navigation (2023/2024)
% Assignment #2: Exercise 3
% Author: Giuseppe Brentino

clearvars; close all; clc;
rng default

addpath('.\kernels\')
addpath('.\sgp4\')
addpath('.\tle\')
addpath('.\mice\src\mice')
addpath('.\mice\lib')
cspice_furnsh('assignment02.tm')
plotStyle;

%% Data
t0 = cspice_str2et('2010-08-12T05:30:00.000');
tf = cspice_str2et('2010-08-12T06:30:00.000');
t_sep = cspice_str2et('August 12 05:27:39.114 UTC 2010');

Mango.TLE = cell(2,1);   
Mango.TLE{1} = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
Mango.TLE{2} = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';
Mango.rsep  = [4622.232026629; 5399.3369588058; -0.0212138165769957];
Mango.vsep  = [0.812221125483763; -0.721512914578826; 7.42665302729053];
Mango.s_sep = [Mango.rsep; Mango.vsep];

Tango.TLE = cell(2,1);   
Tango.TLE{1} = '1 36827U 10028F   10224.22753605  .00278492  00000-0  82287-1 0  9996';
Tango.TLE{2} = '2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55';
Tango.rsep  = [4621.69343340281; 5399.26386352847; -3.09039248714313];
Tango.vsep  = [0.813960847513811; -0.719449862738607; 7.42706066911294];
Tango.s_sep = [Tango.rsep; Tango.vsep];

P_sep = [ +5.6e-7, +3.5e-7, -7.1e-8,    0,        0,        0;
          +3.5e-7, +9.7e-7, +7.6e-8,    0,        0,        0;
          -7.1e-8, +7.6e-8, +8.1e-8,    0,        0,        0;
              0,       0,       0,  +2.8e-11,     0,        0;
              0,       0,       0,      0,    +2.7e-11,     0;
              0,       0,       0,      0,        0,    +9.6e-12 ];

Svalbard.name = 'SVALBARD';
Svalbard.lat = 78.229772;
Svalbard.lon = 15.407786;
Svalbard.alt = 458;
Svalbard.R = diag([0.01^2;0.125^2;0.125^2]); % km deg deg
Svalbard.min_el = 5;

FFRF.R = diag([(1e-5)^2;1;1]); % km deg deg

%% Ex 1

ode_opt = odeset('RelTol',1e-13,'AbsTol',1e-13);
mu = cspice_bodvrd('Earth','GM',1);
J2 = 0.0010826269;
Re = cspice_bodvrd('Earth','RADII',3);
Re = Re(1);

%%% a) Compute visibility window

% propagate orbit
t_span = t0:5:tf;
[~,s1] = ode113(@TBP,[t_sep t0],[Mango.rsep; Mango.vsep],ode_opt,mu);
Mango.s0 = s1(end,:)';
[~,Mango.states_2bp] = ode113(@TBP,t_span,Mango.s0,ode_opt,mu);
Mango.states_2bp = Mango.states_2bp';

% compute first visibility window
[Svalbard.eci,Svalbard.eci2topo] = stationCoordinates (Svalbard.name,t_span);

Svalbard.Mango_coord_2bp = scLocalCoordinates(Svalbard,Mango.states_2bp);

[Svalbard.visibility_time,~] = visibility(Svalbard.Mango_coord_2bp,Svalbard.min_el,t_span);

%%% b) Simulate measurements
[Mango.states_sgp4,Mango.gs_measurements,Svalbard] = simulateGSMeasurements(Mango.TLE,Svalbard);

%%% c) UKF

% Initialization

alpha = 0.1;
beta = 2;
settings_ut.mu = mu;
settings_ut.n = length(Mango.s0);
lambda = alpha^2*settings_ut.n-settings_ut.n;
settings_ut.J2 = J2;
settings_ut.Re = Re;
settings_ut.ode_opt = ode_opt;
settings_ut.lambda = lambda;
settings_ut.fun = 'Absolute';
n = settings_ut.n;

% Define weights
weights.mean = zeros(2*n+1,1);
weights.cov = zeros(2*n+1,1);

weights.mean(1) = lambda/(n+lambda);
weights.cov(1) = weights.mean(1) + (1-alpha^2+beta);

weights.mean(2:end) = 1/(2*(n+lambda)) * ones(2*n,1);
weights.cov(2:end) = weights.mean(2:end);

Mango.gs_measurements = [[t0,0,0,0]; Mango.gs_measurements];  

len = length(Mango.gs_measurements);

% initialize variables
P_abs = zeros(6,6,len);
x_abs = zeros(6,len);
sigma_points = zeros(6,13,len);

% compute initial values
P_sep = 1e4*P_sep;
[sigma_points(:,:,1),x_abs(:,1),P_abs(:,:,1)] = UKFPropagate(t_sep,t0,Mango.s_sep,P_sep,...
    settings_ut,weights);

% initialize performance parameters
sigma_r = zeros(len-1,1);
sigma_v = zeros(len-1,1);
% Sequential filter 

for i=2:len
    Svalbard.eci = Svalbard.eci_vis(:,i-1);
    Svalbard.eci2topo = Svalbard.eci2topo_vis(:,:,i-1);
    [x_abs(:,i),P_abs(:,:,i),sigma_points(:,:,i)] =UKF(x_abs(:,i-1),P_abs(:,:,i-1),...
        settings_ut,Mango.gs_measurements(i-1:i,:),weights,Svalbard);

    sigma_r(i-1) = 3*sqrt(trace(P_abs(1:3,1:3,i)));
    sigma_v(i-1) = 3*sqrt(trace(P_abs(4:6,4:6,i)));
end

% Plot
time_plot = datetime(cspice_timout(Mango.gs_measurements(2:end,1)',...
    'YYYY-MM-DD HR:MN:SC.###'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

error_r = vecnorm(Mango.states_sgp4(1:3,:)-x_abs(1:3,2:end),2,1);
figure()
semilogy(time_plot,error_r)
hold on
semilogy(time_plot,sigma_r)
grid on
xlabel('Time')
ylabel('Range [km]')
legend('error','3$\sigma$')

error_v = vecnorm(Mango.states_sgp4(4:6,:)-x_abs(4:6,2:end),2,1);
figure()
semilogy(time_plot,error_v)
hold on
semilogy(time_plot,sigma_v)
grid on
ylabel('Speed [km/s]')
xlabel('Time')
legend('error','3$\sigma$')

%% Ex 2

%%% a) Compute Tango relative states in LVLH frame

[r_eci, v_eci,~] = TLE2Cartesian(Mango.TLE, t0);
Mango.s0_eci = [r_eci;v_eci];
Mango.mean_motion = sqrt(mu/norm(r_eci)^3);
R_eci2LVLH_0 = eci2Lvlh(Mango.mean_motion,Mango.s0_eci);


[r_eci, v_eci, ~] = TLE2Cartesian(Tango.TLE, t0);
Tango.s0_eci = [r_eci;v_eci];

relative_s0 = R_eci2LVLH_0*(Tango.s0_eci - Mango.s0_eci); % [km],[km/s]

%%% b) Simulate measurements
[Tango.states_window, Tango.sat_measurements_window] = simulateFFRFMeasurements...
    (Mango.mean_motion,t_span,relative_s0,FFRF);

%%% c) Estimate states via UKF

t1 = Mango.gs_measurements(end,1)+5;
t2 = t1+20*60;
t_span = t1:5:t2;
index_i = find(t1==Tango.sat_measurements_window(:,1),1,"first");
index_f = find(t2==Tango.sat_measurements_window(:,1),1,"first");

Tango.states = Tango.states_window(:,index_i:index_f);
Tango.sat_measurements = Tango.sat_measurements_window(index_i:index_f,:);

settings_ut.fun = 'Relative';
settings_ut.mean_motion = Mango.mean_motion;

Tango.sat_measurements = [[t_span(1),0,0,0];Tango.sat_measurements];
len = size(Tango.sat_measurements,1);

P_rel = zeros(6,6,len);
x_rel = zeros(6,len);
sigma_points_rel = zeros(6,13,len);
sigma_r = zeros(len-1,1);
sigma_v = zeros(len-1,1);

P_rel(:,:,1) = diag([0.01, 0.01, 0.1, 0.0001, 0.0001, 0.001]);
x_rel(:,1) = Tango.states(:,1);

for i=2:len
    [x_rel(:,i),P_rel(:,:,i),sigma_points_rel(:,:,i)] =UKF(x_rel(:,i-1),P_rel(:,:,i-1),...
        settings_ut,Tango.sat_measurements(i-1:i,:),weights,FFRF);

    sigma_r(i-1) = 3*sqrt(trace(P_rel(1:3,1:3,i)));
    sigma_v(i-1) = 3*sqrt(trace(P_rel(4:6,4:6,i)));
end

% Plot
time_plot = datetime(cspice_timout(Tango.sat_measurements(2:end,1)',...
    'YYYY-MM-DD HR:MN:SC.###'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
error_r = vecnorm(Tango.states(1:3,:)-x_rel(1:3,2:end),2,1);
figure()
semilogy(time_plot,error_r)
hold on
semilogy(time_plot,sigma_r)
grid on
ylabel('Distance [m]')
xlabel('Time')
legend('error','3$\sigma$')

error_v = vecnorm(Tango.states(4:6,:)-x_rel(4:6,2:end),2,1);
figure()
semilogy(time_plot,error_v)
hold on
semilogy(time_plot,sigma_v)
grid on
ylabel('Speed [m/s]')
xlabel('Time')
legend('error','3$\sigma$')

%% EX 3


time_span = [t_span(1)-5, t_span];
len = length(time_span);
P_prop = zeros(6,6,len);
P_prop(:,:,1) = P_abs(:,:,end);

s_prop = zeros(6,len);
s_prop(:,1) = x_abs(:,1);

P_rot = zeros(6,6,len-1);
Tango.P_abs = zeros(6,6,len-1);
Tango.sigma_r = zeros(len-1,1);
Tango.sigma_v = zeros(len-1,1);

settings_ut.fun = 'Absolute';

for i = 2:len
    %%% a) Propagate Mango covariance
    [~,s_prop(:,i),P_prop(:,:,i)] = UKFPropagate(time_span(i-1),...
        time_span(i),s_prop(:,i-1),P_prop(:,:,i-1),settings_ut,weights);

    %%% b) Rotate Covariance in LVLH frame
    R = eci2Lvlh(Mango.mean_motion,s_prop(:,i));
    P_rot(:,:,i) = R'*P_rel(:,:,i)*R;

    %%% c) Tango absolute covariance
    Tango.P_abs(:,:,i-1) = P_prop(:,:,i-1) + P_rot(:,:,i);
    Tango.sigma_r(i-1) = 3*sqrt(trace(Tango.P_abs(1:3,1:3,i-1))); 
    Tango.sigma_v(i-1) = 3*sqrt(trace(Tango.P_abs(4:6,4:6,i-1)));
end

% Plot
time_plot = datetime(cspice_timout(time_span(2:end),...
    'YYYY-MM-DD HR:MN:SC.###'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
figure()
plot(time_plot,Tango.sigma_r)
ylabel('Distance [km]')
xlabel('Time')
figure()
semilogy(time_plot,Tango.sigma_v)
ylabel('Speed [km/s]')
xlabel('Time')
%% functions

function plotStyle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set figure properties for better looking plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpreter:
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
% lines:
set(0,'defaultLineLineWidth', 1.5);
set(0,'defaultLineMarkerSize',6) ;
set(0,'defaultLineMarkerEdgeColor','k')
set(0,'defaultLineMarkerFaceColor','auto')
% legend:
set(0, 'defaultLegendLocation','southoutside');
set(0, 'defaultLegendOrientation','horizontal');
set(0, 'defaultLegendFontSize',12);
% axes:
set(0,'defaultAxesFontSize',16);
end

function [dxx] = TBP(~,xx,mu)
% TBP - Computes the equations of motion for the Two-Body Problem (TBP)
%       and optionally the state transition matrix (STM) derivative.
%
%   Inputs:
%       - ~: Time variable (not used in this function).
%       - xx: State vector. If the STM is included, xx has 42 elements;
%             otherwise, it has 6 elements.
%       - mu: Gravitational parameter.
%
%   Output:
%       - dxx: Derivative of the state vector. If the STM is included, dxx
%              has 42 elements; otherwise, it has 6 elements.
%

x = xx(1);
y = xx(2);
z = xx(3);
r_norm = norm([x;y;z]);

% Equation of motion
dxx(1:3) = xx(4:6);
dxx(4:6) = -mu/r_norm^3 * [x;y;z];

dxx = dxx';
    
end

function [dxx] = PTBP(t,xx,mu,J2,Re)
% PTBP - Perturbed Two-Body Problem equations of motion with J2 perturbation
%
% Inputs:
%   t   - Time. 
%   xx  - State vector.
%   mu  - Gravitational parameter of the central body.
%   J2  - J2 perturbation coefficient.
%   Re  - Equatorial radius of the central body.
%
% Output:
%   dxx - Derivative of the state vector 

x = xx(1);
y = xx(2);
z = xx(3);
r = [x; y; z];
r_norm = norm(r);

% J2 Perturbation
rotm = cspice_pxform('J2000', 'ITRF93', t);
r_ecef = rotm * r;
a_j2_ecef = 3/2 * J2 * mu * Re^2 / r_norm^5 * r_ecef .* (-[1; 1; 3] + 5 * r_ecef(3)^2 / r_norm^2);
a_j2 = rotm' * a_j2_ecef;

% Equations of motion
dxx(1:3) = xx(4:6);
dxx(4:6) = -mu / r_norm^3 * [x; y; z] + a_j2;

dxx = dxx';

end

function [xx] = CW(t,x0,n)
% CW - Compute the state of a satellite under Clohessy-Wiltshire equations.
%
% Inputs:
%   t  - time between initial and final state
%   x0 - Initial state vector [r; v]
%   n  - Mean motion of the chief satellite
%
% Outputs:
%   xx - State vector [r; v] at time t

% Matrix A components
A_rr = [   4-3*cos(n*t)    0     0;
         6*(sin(n*t)-n*t)  1     0;
                 0         0  cos(n*t) ];

A_vr = [    3*n*sin(n*t)   0       0;
         -6*n*(1-cos(n*t)) 0       0;
                  0        0  -n*sin(n*t) ];

A_rv = [   1/n*sin(n*t)     2/n*(1-cos(n*t))       0;
         -2/n*(1-cos(n*t))  4*sin(n*t)/n-3*t       0 ;
                 0                    0       1/n*sin(n*t) ];

A_vv = [   cos(n*t)    2*sin(n*t)      0;
         -2*sin(n*t)  4*cos(n*t)-3     0;
               0            0       cos(n*t) ];

% Assemble matrix
A = [ A_rr A_rv;
      A_vr A_vv ];

% Compute the state
xx = (A*x0)';

end

function [station_eci, station_eci2topo] = stationCoordinates(station_name, t_span)
% stationCoordinates - Compute station coordinates in ECI and
% transformation to topocentric frame
%
% Inputs:
%   station_name  - Station name (e.g., 'STATION1').
%   t_span        - Time span for which coordinates are computed
%
% Outputs:
%   station_eci       - Station coordinates in ECI frame
%   station_eci2topo  - Transformation matrix from ECI to topocentric frame

% Define topocentric frame name based on station name
station_topo_frame = [station_name, '_TOPO'];

% Transformation from ECI to topocentric frame
station_eci2topo = cspice_sxform('J2000', station_topo_frame, t_span);

% Compute station position in ECI
station_eci = cspice_spkezr(station_name, t_span, 'J2000', 'NONE', 'EARTH');

end

function [sat_coord,varargout] = scLocalCoordinates(station,states,varargin)
% scLocalCoordinates - Compute topocentric coordinates of a satellite with
% respect to a ground station.
%
% Inputs:
%   station - Structure containing information about the ground station
%   states  - Satellite states in ECI coordinates (position and velocity)
%   varargin - Additional optional arguments (visibility time window)
%
% Outputs:
%   sat_coord - Matrix containing range, azimuth, and elevation
%   varargout - Optional output (visibility information, measures and
%               updated station structure)

% Compute station-satellite vector in ECI
sat_eci = states - station.eci;
n = size(sat_eci, 2);
sat_topo = zeros(6, n);

% Convert state into topocentric frame
for i = 1:n
    sat_topo(:, i) = station.eci2topo(:, :, i) * sat_eci(:, i);
end

% Compute range, azimuth, and elevation using cspice_xfmsta
sat_coord = cspice_xfmsta(sat_topo, 'RECTANGULAR', 'LATITUDINAL', 'EARTH');
sat_Az = wrapTo360(sat_coord(2, :) * cspice_dpr);
sat_El = sat_coord(3, :) * cspice_dpr;
sat_range = sat_coord(1, :);
sat_coord = [sat_range; sat_Az; sat_El];

% Simulating measurements and handling visibility
if nargout == 4

    if nargin == 3

        % Simulate measurements with Gaussian noise
        measures = mvnrnd(sat_coord(1:3, :)', station.R);
        visibility = measures(:, 3) >= station.min_el;

        % Update station information based on visibility
        station.eci_vis = station.eci(:, visibility);
        station.eci2topo_vis = station.eci2topo(:, :, visibility);

    measures = measures(visibility,:);
        vis_time = varargin{1}(visibility)';

        % Output visibility information, measures and updated station structure
        varargout{1} = [vis_time, measures];
        varargout{2} = station;
        varargout{3} = visibility;

    else
        error('Visibility window required to compute measurements');
    end

end

end

function [sat_coord, varargout] = scRelCoordinates(states, varargin)
% scRelCoordinates - Convert Cartesian coordinates to spherical coordinates.
%
% Inputs:
%   states   - Matrix of Cartesian coordinates (position and velocity)
%   varargin - Optional input (measurement noise covariance matrix)
%
% Outputs:
%   sat_coord - Matrix containing range, azimuth, and elevation
%   varargout - Optional output (noisy spherical coordinates)

% Convert from cartesian to spherical coordinates
sat_coord = cspice_xfmsta(states, 'RECTANGULAR', 'LATITUDINAL', 'Earth');
sat_Az = wrapTo360(sat_coord(2, :) * cspice_dpr);
sat_El = sat_coord(3, :) * cspice_dpr;
sat_range = sat_coord(1, :);
sat_coord = [sat_range; sat_Az; sat_El];

% Add noise
if nargout == 2
    R = varargin{1};

    % Simulate measurements with Gaussian noise
    varargout{1} = mvnrnd(sat_coord', R);
end

end

function [visibility_time, coord] = visibility(sat_coord, min_el, t_span)
% visibility - Extract visibility time and coordinates based on minimum elevation.
%
% Inputs:
%   sat_coord   - Matrix containing range, azimuth, elevation
%   min_el      - Minimum elevation for visibility
%   t_span      - Time span
%
% Outputs:
%   visibility_time - Extracted visibility times
%   coord           - Corresponding coordinates

% Extract visibility condition based on minimum elevation
visibility_condition = sat_coord(3, :) >= min_el;

% Extract visibility time and coordinates
visibility_time = t_span(visibility_condition);
coord = sat_coord(:, visibility_condition);

end

function [r_eci, v_eci, initial_epoch] = TLE2Cartesian(TLE, varargin)
% TLE2Cartesian - Convert Two-Line Element (TLE) to Cartesian state vectors
%
% Inputs:
%   TLE          - Two-Line Element data in a cell array, woth one element
%                  per line
%   varargin     - Optional input for time span in ET (in seconds)
%
% Outputs:
%   r_eci        - ECI position vectors (3xN)
%   v_eci        - ECI velocity vectors (3xN)
%   initial_epoch - Initial epoch of the TLE data in ET

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% Convert to rv
initial_data = twoline2rv(TLE{1}, TLE{2}, typerun, 'e', opsmode, whichconst);

% Get TLE epoch
[year, mon, day, hr, min, sec] = invjday(initial_data.jdsatepoch, initial_data.jdsatepochf);
initial_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year, mon, day, hr, min, sec]);
initial_epoch = cspice_str2et(initial_epoch_str);

% Compute s/c states during visibility windows
if nargin == 1
    t_span = initial_epoch;
    t_span_TLE = 0;
else
    t_span = varargin{1};
    t_span_TLE = (t_span - initial_epoch) / 60;  % Convert time span to minutes
end

n = length(t_span_TLE);
r_teme = zeros(3, n);
v_teme = zeros(3, n);
ateme = [0; 0; 0];

% Nutation correction
% delta-Psi and delta-Epsilon corrections on the observation day,
% from https://celestrak.org/SpaceData/EOP-All.csv
dpsi = -0.073296 * pi / (180 * 3600); %[rad]
deps = -0.009373 * pi / (180 * 3600); %[rad]

r_eci = zeros(3, length(t_span_TLE));
v_eci = zeros(3, length(t_span_TLE));

for i = 1:length(t_span_TLE)

    % Compute states in TEME reference frame
    [~, r_teme(:, i), v_teme(:, i)] = sgp4(initial_data, t_span_TLE(i));

    % Precession correction: centuries from TDT 2000 January 1 00:00:00.000
    ttt = cspice_unitim(t_span(i), 'ET', 'TDT') / cspice_jyear() / 100;

    % Rotation from TEME to ECI
    [r_eci(:, i), v_eci(:, i), ~] = teme2eci(r_teme(:, i), v_teme(:, i), ateme, ttt, dpsi, deps);
end

end

function [states, measure,station] = simulateGSMeasurements(TLE, station)
% SimulateGSMeasurements - Simulate measurements of a satellite position from
% a ground station
%
% Inputs:
%   TLE      - Two-Line Element data in a cell array, woth one element per
%              line
%   station  - Structure containing ground station information
%
% Outputs:
%   coord    - Satellite coordinates in the topocentric frame (range, 
%              azimuth, elevation)
%   measure  - Simulated measurements (range, azimuth, elevation) with
%              added noise

% Simulate measurements
[station.eci, station.eci2topo] = stationCoordinates(station.name, station.visibility_time);
[r_eci, v_eci, ~] = TLE2Cartesian(TLE, station.visibility_time);
states = [r_eci;v_eci];
[~, measure,station,vis_flag] = scLocalCoordinates(station, states,station.visibility_time);
states = states(:,vis_flag);

end

function [states, measures] = simulateFFRFMeasurements(n, t_span, x0, FFRF)
% simulateFFRFMeasurements - Simulate relative states and measurements in 
% the FFRF sensor frame (assumed coincident with the satellite frame).
%
% Inputs:
%   n      - Mean motion of the chief satellite
%   t_span - Time span
%   x0     - Initial state vector
%   FFRF   - Structure containing FFRF-related information
%
% Outputs:
%   states   - Matrix containing simulated states over time
%   measures - Matrix containing simulated measurements (with noise) and
%              corresponding timestamps

% Initialize states
states = zeros(6, length(t_span));
states(:, 1) = x0;

% Simulate states using Clohessy-Wiltshire equations
for i = 2:length(t_span)
    tof = t_span(i) - t_span(i-1);
    states(:, i) = CW(tof, states(:, i - 1), n);
end

% Simulate measurements in the FFRF frame
[~, measures] = scRelCoordinates(states, FFRF.R);
measures = [t_span', measures];

end

function R = eci2Lvlh(n,s)
% eci2Lvlh - Compute the transformation matrix from ECI to LVLH frame.
%
% Inputs:
%   n - Orbiting body's mean motion [rad/s]
%   s - State vector in ECI frame [position; velocity]
%
% Output:
%   R - Transformation matrix from ECI to LVLH frame

% Extract position and velocity vectors from the state
r_eci = s(1:3);
v_eci = s(4:6);

% Calculate the unit vectors in the LVLH frame
i_v = r_eci / norm(r_eci);
k_v = cross(r_eci, v_eci) / norm(cross(r_eci, v_eci));
j_v = cross(k_v, i_v);

% Construct the Ssew-symmetric matrix
S = [0  n  0;
    -n  0  0;
     0  0  0];

% Assemble the transformation matrix from ECI to LVLH
R = [ [i_v'; j_v'; k_v']             zeros(3);
      S * [i_v'; j_v'; k_v']    [i_v'; j_v'; k_v'] ];

end

function [x_hat, P, sigma_points] = UKF(x_hat, P, settings_ut, measure, weights, station)
% UKF - Unscented Kalman Filter for state estimation. This function serves
%       as a wrapper for the UKF process, including the prediction
%       (propagation) step (UKFPropagate) and the correction step
%       (UKFCorrect).
%
% Inputs:
%   x_hat - State estimate at the previous time step
%   P - Covariance matrix of the state estimate at the previous time step
%   settings_ut - Settings for the Unscented Transform (UT)
%   measure - Measurement data at the current time step
%   weights - Weights used in the UKF calculations
%   station - Information about the tracking station
%
% Outputs:
%   x_hat - Updated state estimate
%   P - Updated covariance matrix of the state estimate
%   sigma_points - Sigma points used in the filtering process

% Extract the initial and final time from the measurement data
t0 = measure(1, 1);
tf = measure(2, 1);

% Propagate sigma points through the process model
[sigma_points, x_hat, P] = UKFPropagate(t0, tf, x_hat, P, settings_ut, weights);

% Correct the state estimate and covariance using the measurement
[x_hat, P] = UKFCorrect(x_hat, P, sigma_points, measure, station, weights, settings_ut);

end

function [sigma_points, x_hat, P] = UKFPropagate(t0, tf, x_hat, P, settings_ut, weights)
% UKFPropagate - Propagate sigma points and covariance matrix with 
%                Unscented Transform.
%
% Inputs:
%   t0 - Initial time
%   tf - Final time
%   x_hat - State estimate at the previous time step
%   P - Covariance matrix of the state estimate at the previous time step
%   settings_ut - Settings for the Unscented Transform 
%   weights - Weights used in the UKF calculations
%
% Outputs:
%   sigma_points - Propagated sigma points
%   x_hat - Updated state estimate after propagation
%   P - Updated covariance matrix of the state estimate after propagation

n = settings_ut.n; % Dimension of the state vector
ode_opt = settings_ut.ode_opt; % ODE solver options

lambda = settings_ut.lambda; % Scaling factor

% Compute sigma points at t_(k-1)
mat = sqrtm((n + lambda) * P);
sigma_points(:, 1, 1) = x_hat;
for i = 1:n
    sigma_points(:, i + 1, 1) = x_hat + mat(:, i);
    sigma_points(:, i + 1 + n, 1) = x_hat - mat(:, i);
end

% Propagate sigma points to t_k
if t0 < tf
    for i = 1:(2 * n + 1)
        switch settings_ut.fun

            % Compute absolute state using the perturbed two body problem
            case 'Absolute'            
                mu = settings_ut.mu;
                J2 = settings_ut.J2;
                Re = settings_ut.Re;
                [~, x] = ode113(@PTBP, [t0, tf], sigma_points(:, i), ode_opt, mu, J2, Re);

            % Compute relative state using the Clohessy-Wiltshire equations
            case 'Relative'
                mean_motion = settings_ut.mean_motion;
                x = CW(tf - t0, sigma_points(:, i), mean_motion);
        end
        sigma_points(:, i) = x(end, :)';
    end

    % Compute state vector
    x_hat = zeros(6, 1);
    for i = 1:size(sigma_points, 2)
        x_hat = x_hat + weights.mean(i) * sigma_points(:, i);
    end

    % Compute Covariance matrix
    P = zeros(6, 6);
    for i = 1:size(sigma_points, 2)
        P = P + weights.cov(i) * (sigma_points(:, i) - x_hat) * (sigma_points(:, i) - x_hat)';
    end
end

end

function [x_hat, P] = UKFCorrect(x_hat, P, sigma_points, measure, station, weights, settings_ut)
% UKFCorrect - Correct the state estimate based on measurements.
%
% Inputs:
%   x_hat - State estimate after propagation
%   P - Covariance matrix of the state estimate after propagation
%   sigma_points - Propagated sigma points
%   measure - Measurement vector
%   station - Tracking station information
%   weights - Weights used in the UKF calculations
%   settings_ut - Settings for the Unscented Transform
%
% Outputs:
%   x_hat - Corrected state estimate
%   P - Updated covariance matrix of the state estimate

n = settings_ut.n; % Dimension of the state vector

% Simulate measurements using propagated sigma points
gamma = zeros(3, 2 * n + 1);
for i = 1:size(gamma, 2)

    switch settings_ut.fun
        case 'Absolute'
            gamma(:, i) = scLocalCoordinates(station, sigma_points(:, i));
        case 'Relative'
            gamma(:, i) = scRelCoordinates(sigma_points(:, i));
    end

end

% Compute the predicted measurement
y_hat = zeros(3, 1);
for i = 1:size(gamma, 2)
    y_hat = y_hat + weights.mean(i) * gamma(:, i);
end

% Covariance matrices

Pyy = station.R;
meas_diff = zeros(3, size(gamma, 2));
for i = 1:size(gamma, 2)
    meas_diff(:, i) = [gamma(1, i)' - y_hat(1); ...
        180/pi * angdiff(y_hat(2:3) * pi/180, gamma(2:3, i) * pi/180)];
    Pyy = Pyy + weights.cov(i) * meas_diff(:, i) * meas_diff(:, i)';
end

Pxy = zeros(6, 3);
for i = 1:size(gamma, 2)
    Pxy = Pxy + weights.cov(i) * (sigma_points(:, i) - x_hat) * meas_diff(:, i)';
end

% Perform correction if Pyy is invertible
if cond(Pyy) < 1e10

    % Kalman gain
    K = Pxy / Pyy;

    % Correction of the state
    e = [measure(2, 2)' - y_hat(1); ...
        180/pi * angdiff(y_hat(2:3) * pi/180, measure(2, 3:4)' * pi/180)];
    x_hat = x_hat + K * e;

    % Propagation of covariance matrix
    P = P - K * Pyy * K';
end

end

