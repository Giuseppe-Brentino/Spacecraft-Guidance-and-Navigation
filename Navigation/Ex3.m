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
Mango.r0  = [4622.232026629; 5399.3369588058; -0.0212138165769957];
Mango.v0  = [0.812221125483763; -0.721512914578826; 7.42665302729053];

Mango.s_sep = [Mango.r0; Mango.v0];

P0 = [ +5.6e-7, +3.5e-7, -7.1e-8,    0,        0,        0;
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
t_sep= cspice_str2et('August 12 05:27:39.114 UTC 2010');

ode_opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
mu = cspice_bodvrd('Earth','GM',1);
J2 = 0.0010826269;
Re = cspice_bodvrd('Earth','RADII',3);
Re = Re(1);

%%% a)

% propagate orbit
t_span = t0:5:tf;
[~,s1] = ode113(@PTBP,[t_sep t0],[Mango.r0; Mango.v0],ode_opt,mu,J2,Re);
Mango.s0 = s1(end,:)';
[~,Mango.states_2bp] = ode113(@PTBP,t_span,Mango.s0,ode_opt,mu,J2,Re);
Mango.states_2bp = Mango.states_2bp';

% compute first visibility window
[Svalbard.eci,Svalbard.eci2topo] = stationCoordinates (Svalbard,t_span);

Svalbard.Mango_coord_2bp = scLocalCoordinates(Svalbard,Mango.states_2bp);

[visibility_window,~] = visibility(Svalbard.Mango_coord_2bp,t_span,Svalbard.min_el);

%%% b) Simulate measurements
[Mango.states_sgp4,Mango.gs_measurements,Svalbard] = simulateGSMeasurements(Mango.TLE,Svalbard,visibility_window);

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
mat = sqrtm((n+lambda)*P0);
    sigma_points(:,1,1) = Mango.s_sep;
    for i = 1:n
        sigma_points(:,i+1,1) = Mango.s_sep + mat(:,i);
        sigma_points(:,i+1+n,1) = Mango.s_sep - mat(:,i);
    end

% propagate from separation time to t0
[sigma_points(:,:,1), P(:,:,1), x(:,1)] = UT(t_sep,t0,sigma_points(:,:,1),...
    settings_ut,weights);
P(:,:,1) = P(:,:,1)*1e4;



% initialize performance parameters
sigma_r = zeros(len-1,1);
sigma_v = zeros(len-1,1);
% Sequential filter 

for i=2:len
    Svalbard.eci = Svalbard.eci_vis(:,i-1);
    Svalbard.eci2topo = Svalbard.eci2topo_vis(:,:,i-1);
    [x(:,i),P(:,:,i),sigma_points(:,:,i)] =UKF(sigma_points(:,:,i-1),...
        settings_ut,Mango.gs_measurements(i-1:i,:),weights,Svalbard);

    sigma_r(i-1) = 3*sqrt(trace(P(1:3,1:3,i)));
    sigma_v(i-1) = 3*sqrt(trace(P(4:6,4:6,i)));
end

% Plot
figure()
hold on 
grid on
plot(Mango.gs_measurements(2:end,1),Mango.states_sgp4(3,:))
plot(Mango.gs_measurements(2:end,1),x(3,2:end),'--')
legend('SGP4','UKF')

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

function [dxx] = PTBP(t,xx,mu,J2,Re)

x = xx(1);
y = xx(2);
z = xx(3);
r = [x;y;z];
r_norm = norm(r);

%J2
rotm= cspice_pxform('J2000', 'ITRF93', t);
r_ecef = rotm*r;

a_j2_ecef = 3/2 * J2 * mu * Re^2/r_norm^5 * r_ecef.*( -[1;1;3] + 5*r_ecef(3)^2/r_norm^2); 

a_j2 = rotm'*a_j2_ecef;

% RHS
dxx(1:3) = xx(4:6);
dxx(4:6) = -mu/r_norm^3 * [x;y;z] + a_j2;

dxx = dxx';

end

function [station_eci,station_eci2topo] = stationCoordinates (station,t_span)

% Define station name
station.topo_frame = [station.name, '_TOPO'];

% Transformation from ECI to topocentric frame
station_eci2topo = cspice_sxform('J2000', station.topo_frame, t_span);

% Compute station position in ECI
station_eci = cspice_spkezr(station.name, t_span, 'J2000', 'NONE', 'EARTH');

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
if nargout == 3 
    if nargin == 3
    measures = mvnrnd(sat_coord(1:3,:)',station.R);
    visibility = measures(:,3)>=station.min_el;
    station.eci_vis = station.eci(:,visibility);
    station.eci2topo_vis = station.eci2topo(:,:,visibility);
    measures = measures(visibility,:);
    vis_time = varargin{1}(visibility)';
    varargout{1} = [vis_time, measures];
   
    varargout{2} = station;
    else
        error('Visibility window required to compure measurements')
    end
end

end

function [visibility_time, coord] = visibility(sat_coord,t_span,min_el)
    
    visibility_condition = sat_coord(3,:)>=min_el;
    index_i = find(visibility_condition==true,1,'first');
    index_f = index_i + find(visibility_condition(index_i:end)==false,1,'first');
    visibility_time = t_span(index_i:index_f-1);
    coord = sat_coord(:,index_i:index_f-1);

end

function [states,measure,station] = simulateGSMeasurements(TLE,station,t_span)

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% convert to rv
sat.initial_data = twoline2rv(TLE{1}, TLE{2}, typerun,'e', opsmode, whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(sat.initial_data.jdsatepoch, sat.initial_data.jdsatepochf);
sat.initial_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat.initial_epoch = cspice_str2et(sat.initial_epoch_str);

%%%% Compute s/c states during visibility windows

sat.time_TLE = (t_span - sat.initial_epoch)./60;

n = length(sat.time_TLE);
sat.r_teme = zeros(3,n);
sat.v_teme = zeros(3,n);
ateme = [0;0;0];

% Nutation correction
% delta-Psi and delta-Epsilon corrections on the observation day, from https://celestrak.org/SpaceData/EOP-All.csv
dpsi = -0.073296 * pi / (180*3600); %[rad]
deps = -0.009373 * pi / (180*3600); %[rad]


for i = 1:length(sat.time_TLE)
    % Compute states in TEME reference frame 
[~,sat.r_teme(:,i),sat.v_teme(:,i)] = sgp4(sat.initial_data,sat.time_TLE(i));

% Precession correction: centuries from TDT 2000 January 1 00:00:00.000
ttt = cspice_unitim(t_span(i), 'ET', 'TDT')/cspice_jyear()/100;

% rotation from TEME to ECI
[sat.r_eci(:,i),sat.v_eci(:,i),~] = teme2eci(sat.r_teme(:,i),sat.v_teme(:,i), ateme, ttt, dpsi,deps);
end

% simulate measurements
states = [sat.r_eci;sat.v_eci];
[station.eci,station.eci2topo] = stationCoordinates (station,t_span);
[~,measure,station] = scLocalCoordinates(station,states,t_span);
states = states(:,ismember(t_span,measure(:,1)'));
end

function [x_hat,P,sigma_points] = UKF(sigma_points,settings_ut,measure,weights,station)

n = settings_ut.n;

t0 = measure(1,1);
tf = measure(2,1);

%%% propagation step
  
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
        K = Pxy/(Pyy);

        % Correction of the state
        e = [measure(2,2)'-y_hat(1);...
            180/pi*angdiff(measure(2,3:4)'*pi/180, y_hat(2:3)*pi/180) ];
        x_hat = x_hat + K*e;

        % Propagation of covariance matrix
        P = P*K*Pyy*K';
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