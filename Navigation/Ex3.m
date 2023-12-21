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

Mango.TLE = cell(2,1);   
Mango.TLE{1} = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
Mango.TLE{2} = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';
Mango.rsep  = [4622.232026629; 5399.3369588058; -0.0212138165769957];
Mango.vsep  = [0.812221125483763; -0.721512914578826; 7.42665302729053];

Mango.s_sep = [Mango.rsep; Mango.vsep];

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

%% Ex 1

t0 = cspice_str2et('2010-08-12T05:30:00.000');
tf = cspice_str2et('2010-08-12T06:30:00.000');
t_sep = cspice_str2et('August 12 05:27:39.114 UTC 2010');

ode_opt = odeset('RelTol',1e-12,'AbsTol',[ones(3,1)*1e-9; ones(3,1)*1e-12]);
mu = cspice_bodvrd('Earth','GM',1);
J2 = 0.0010826269;
Re = cspice_bodvrd('Earth','RADII',3);
Re = Re(1);

%%% a)

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
P = zeros(6,6,len);
x = zeros(6,len);
sigma_points = zeros(6,13,len);

% compute initial values
P_sep = 1e4*P_sep;
mat = sqrtm((n+lambda)*P_sep);
    sigma_points(:,1,1) = Mango.s_sep;
    for i = 1:n
        sigma_points(:,i+1,1) = Mango.s_sep + mat(:,i);
        sigma_points(:,i+1+n,1) = Mango.s_sep - mat(:,i);
    end

% propagate from separation time to t0
[sigma_points(:,:,1), P(:,:,1), x(:,1)] = UT(t_sep,t0,sigma_points(:,:,1),...
    settings_ut,weights);

% initialize performance parameters
sigma_r = zeros(len-1,1);
sigma_v = zeros(len-1,1);
% Sequential filter 

for i=2:len
    Svalbard.eci = Svalbard.eci_vis(:,i-1);
    Svalbard.eci2topo = Svalbard.eci2topo_vis(:,:,i-1);
    [x(:,i),P(:,:,i),sigma_points(:,:,i)] =UKF(x(:,i-1),P(:,:,i-1),...
        settings_ut,Mango.gs_measurements(i-1:i,:),weights,Svalbard);

    sigma_r(i-1) = 3*sqrt(trace(P(1:3,1:3,i)));
    sigma_v(i-1) = 3*sqrt(trace(P(4:6,4:6,i)));
end

% Plot
error_r = vecnorm(Mango.states_sgp4(1:3,:)-x(1:3,2:end),2,1);
figure()
semilogy(Mango.gs_measurements(2:end,1),error_r)
hold on
semilogy(Mango.gs_measurements(2:end,1),sigma_r)
grid on
title('Position')
legend('error','3$\sigma$')

error_v = vecnorm(Mango.states_sgp4(4:6,:)-x(4:6,2:end),2,1);
figure()
semilogy(Mango.gs_measurements(2:end,1),error_v)
hold on
semilogy(Mango.gs_measurements(2:end,1),sigma_v)
grid on
title('Velocity')
legend('error','3$\sigma$')
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

% Compute station-satellite vector in ECI
sat_eci = states - station.eci;
n = size(sat_eci,2);
sat_topo = zeros(6,n);
% Convert state into topocentric frame
for i = 1:n
   sat_topo(:,i) = station.eci2topo(:,:,i)*sat_eci(:,i);
end

% Compute range, azimuth and elevation using cspice_xfmsta
sat_coord = cspice_xfmsta( sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');
sat_Az = wrapTo360( sat_coord(2,:)*cspice_dpr );
sat_El = sat_coord(3,:)*cspice_dpr;
sat_range = sat_coord(1,:);
sat_coord = [sat_range;sat_Az;sat_El];
if nargout == 4 
    if nargin == 3
    measures = mvnrnd(sat_coord(1:3,:)',station.R);
    visibility = measures(:,3)>=station.min_el;
    station.eci_vis = station.eci(:,visibility);
    station.eci2topo_vis = station.eci2topo(:,:,visibility);
    measures = measures(visibility,:);
    vis_time = varargin{1}(visibility)';
    varargout{1} = [vis_time, measures];
    varargout{2} = station;
    varargout{3} = visibility;
    else
        error('Visibility window required to compure measurements')
    end
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
% delta-Psi and delta-Epsilon corrections on the observation day, from https://celestrak.org/SpaceData/EOP-All.csv
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

function [x_hat,P,sigma_points] = UKF(x_hat,P,settings_ut,measure,weights,station)

n = settings_ut.n;
lambda = settings_ut.lambda;

t0 = measure(1,1);
tf = measure(2,1);

%%% propagation step
 mat = sqrtm((n+lambda)*P);
    sigma_points(:,1,1) = x_hat;
    for i = 1:n
        sigma_points(:,i+1,1) = x_hat + mat(:,i);
        sigma_points(:,i+1+n,1) = x_hat - mat(:,i);
    end

[sigma_points, P, x_hat] = UT(t0,tf,sigma_points, settings_ut,weights);

    %%% correction step
    % measurements
    gamma = zeros(3,2*n+1);
    for i = 1:size(gamma,2)
        gamma(:,i) = scLocalCoordinates(station,sigma_points(:,i));
    end
    y_hat = zeros(3,1);
    for i = 1:size(gamma,2)
        y_hat = y_hat + weights.mean(i)*gamma(:,i);
    end
    % covariance matrices
    Pyy = station.R;
    for i = 1:size(gamma,2)
        Pyy = Pyy + weights.cov(i) * (gamma(:,i) - y_hat) * (gamma(:,i) - y_hat)';
    end

    Pxy = zeros(6,3);
    for i = 1:size(gamma,2)
        Pxy = Pxy +  weights.cov(i) * (sigma_points(:,i) -x_hat) * (gamma(:,i) -y_hat)';
    end

    if cond(Pyy) < 1e4 % perform correction if Pyy is invertible

        % kalman gain
        K = Pxy/Pyy;
        % Correction of the state
        e = [measure(2,2)'-y_hat(1);...
            180/pi*angdiff( y_hat(2:3)*pi/180, measure(2,3:4)'*pi/180) ];
        x_hat = x_hat + K*e;
        % Propagation of covariance matrix
        P = P-K*Pyy*K';
    end
end

function [sigma_points, P, x_hat] = UT(t0,tf,sigma_points, settings_ut,weights)

mu = settings_ut.mu;
J2 = settings_ut.J2;
Re = settings_ut.Re;
n = settings_ut.n;
ode_opt = settings_ut.ode_opt;

%sigma points
for i = 1:2*n+1
    [~,x] = ode113(@PTBP,[t0,tf],sigma_points(:,i),ode_opt,mu,J2,Re);
    sigma_points(:,i) = x(end,:)';
end

% state
x_hat = zeros(6,1);
for i = 1:size(sigma_points,2)
    x_hat = x_hat + weights.mean(i)*sigma_points(:,i);
end

% covariance matrix
P = zeros(6,6);
for i = 1:size(sigma_points,2)
    P = P + weights.cov(i) * (sigma_points(:,i) -x_hat) * (sigma_points(:,i) -x_hat)';
end
end
