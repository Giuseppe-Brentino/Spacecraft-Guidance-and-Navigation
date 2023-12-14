% Spacecraft Guidance and Navigation (2023/2024)
% Assignment #2: Exercise 1
% Author: Giuseppe Brentino

clearvars; close all; clc;

addpath('.\kernels\')
addpath('.\sgp4\')
addpath('.\tle\')
addpath('.\mice\src\mice')
addpath('.\mice\lib')
cspice_furnsh('assignment02.tm')
plotStyle;

%% Init data

Kourou.name = 'KOUROU';
Kourou.lat = 5.25144;
Kourou.lon = -52.80466;
Kourou.alt = -14.67;
Kourou.R = diag([0.01^2;0.1^2;0.1^2]); % km deg deg (^2)
Kourou.W_m = inv(sqrtm(Kourou.R));
Kourou.min_el = 10;

Svalbard.name = 'SVALBARD';
Svalbard.lat = 78.229772;
Svalbard.lon = 15.407786;
Svalbard.alt = 458;
Svalbard.R = diag([0.01^2;0.125^2;0.125^2]); % km deg deg
Svalbard.W_m  = inv(sqrtm(Svalbard.R));
Svalbard.min_el = 5;

r0  = [4622.232026629; 5399.3369588058; -0.0212138165769957];
v0  = [0.812221125483763; -0.721512914578826; 7.42665302729053]; 
t0 = cspice_str2et('2010-08-12T05:30:00.000');
tf = cspice_str2et('2010-08-12T11:00:00.000');
dt = 60;
t_span = t0:dt:tf;

J2 = 0.0010826269;
Re = cspice_bodvrd('Earth','RADII',3);
Re = Re(1);
%% Ex 1
mu = cspice_bodvrd('Earth','GM',1);
ode_opt = odeset('RelTol',1e-13,'AbsTol',1e-13);

[~,states] = ode78(@TBP,t_span,[r0;v0;],ode_opt,mu);
states = states';

[Kourou.eci,Kourou.eci2topo] = stationCoordinates (Kourou,t_span);
[Svalbard.eci,Svalbard.eci2topo] = stationCoordinates (Svalbard,t_span);

Kourou.Mango_coord = scLocalCoordinates(Kourou,states);
Svalbard.Mango_coord = scLocalCoordinates(Svalbard,states);

[Kourou.visibility_time,Kourou.Mango_visible_coord] = ...
    visibility(Kourou.Mango_coord,Kourou.min_el,t_span);
[Svalbard.visibility_time,Svalbard.Mango_visible_coord] = ...
    visibility(Svalbard.Mango_coord,Svalbard.min_el,t_span);

% plot visibility coordinates
figure()
% Kourou
subplot(2,2,1)
title('Kourou')
hold on
grid on
plot(t_span,Kourou.Mango_coord(2,:),'DisplayName','Not visible');

vis_time = Kourou.visibility_time;
Az_plot = Kourou.Mango_visible_coord(2,:);
El_plot = Kourou.Mango_visible_coord(3,:);
for i = 2:length(vis_time)
    if vis_time(i) - vis_time(i-1) >= 2*dt
        vis_time =[vis_time(1:i-1) nan vis_time(i:end)];
        Az_plot = [Az_plot(1:i-1) nan Az_plot(i:end)];
        El_plot = [ El_plot(1:i-1) nan  El_plot(i:end)];
    end
end
plot(vis_time,Az_plot,'DisplayName','Visible');
ylabel('Azimuth [deg]')
xlabel('Epoch [s]')
legend
subplot(2,2,3)
hold on
grid on
plot(t_span,Kourou.Mango_coord(3,:),'DisplayName','Not visible');
plot(vis_time,El_plot,'DisplayName','Visible');
ylabel('Elevation [deg]')
xlabel('Epoch [s]')
legend;

%Svalbard
subplot(2,2,2)
title('Svalbard')
hold on
grid on
plot(t_span,Svalbard.Mango_coord(2,:),'DisplayName','Not visible');

vis_time = Svalbard.visibility_time;
Az_plot = Svalbard.Mango_visible_coord(2,:);
El_plot = Svalbard.Mango_visible_coord(3,:);
for i = 2:length(vis_time)
    if vis_time(i) - vis_time(i-1) >= 2*dt
        vis_time =[vis_time(1:i-1) nan vis_time(i:end)];
        Az_plot = [Az_plot(1:i-1) nan Az_plot(i:end)];
        El_plot = [El_plot(1:i-1) nan  El_plot(i:end)];
    end
end

plot(vis_time,Az_plot,'DisplayName','Visible');
ylabel('Azimuth [deg]')
xlabel('Epoch [s]')
legend
subplot(2,2,4)
hold on
grid on
plot(t_span,Svalbard.Mango_coord(3,:),'DisplayName','Not visible');
plot(vis_time,El_plot,'DisplayName','Visible');
ylabel('Elevation [deg]')
xlabel('Epoch [s]')
legend;

%% Ex 2

% Initial TLE
TLE = cell(2,1);   
TLE{1} = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
TLE{2} = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';

[Kourou,Svalbard,Mango,Kourou.Mango_measured_coord,Svalbard.Mango_measured_coord] = ...
    simulateMeasurements(TLE,Kourou,Svalbard);


%% Ex 3

% a)
fun =@(t,x)TBP(t,x,mu);
x0_guess = [Mango.r_eci(:,1); Mango.v_eci(:,1)] ;
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display', 'iter');

data = cell(1,2);
data{1,1} = Kourou;
data{1,2} = Kourou.Mango_measured_coord;
[x0_Kourou,resnorm_Kourou,residual_Kourou,exitflag_Kourou,~,~,jacobian_Kourou] = ...
lsqnonlin(@costFcn, x0_guess, [],[],options,data,mu,fun,Mango.initial_epoch,ode_opt);

% b)
data = cell(2,2);
data{1,1} = Kourou;
data{1,2} = Kourou.Mango_measured_coord;
data{2,1} = Svalbard;
data{2,2} = Svalbard.Mango_measured_coord;
[x0_both,resnorm_both,residual_both,exitflag_both,~,~,jacobian_both] = ...
lsqnonlin(@costFcn, x0_guess, [],[],options,data,mu,fun,Mango.initial_epoch,ode_opt);

% c)
fun =@(t,x)PTBP(t,x,mu,J2,Re);
[x0_J2,resnorm_J2,residual_J2,exitflag_J2,~,~,jacobian_J2] = ...
lsqnonlin(@costFcn, x0_guess, [],[],options,data,mu,fun,Mango.initial_epoch,ode_opt,J2,Re);

sol1 = x0_guess - x0_Kourou
sol2 = x0_guess - x0_both
sol3 = x0_guess - x0_J2

%% Ex 5

% Covariance computation
P_ls_Kourou = resnorm_Kourou/(length(residual_Kourou)-length(x0_Kourou))...
    .*inv(jacobian_Kourou.'*jacobian_Kourou);
P_ls_both = resnorm_both/(length(residual_both)-length(x0_both))...
    .*inv(jacobian_both.'*jacobian_both);
P_ls_J2 = resnorm_J2/(length(residual_J2)-length(x0_J2))...
    .*inv(jacobian_J2.'*jacobian_J2);

% Best combination : J2(better model) + both station (more measurements)

% Initial TLE
TLE = cell(2,1);   
TLE{1} = '1 36827U 10028F   10224.22753605  .00278492  00000-0  52287-1 0  9996';
TLE{2} = '2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55';

[Kourou,Svalbard,Tango,Kourou.Tango_measured_coord,Svalbard.Tango_measured_coord] = ...
    simulateMeasurements(TLE,Kourou,Svalbard);

x0_guess = [Tango.r_eci(:,1); Tango.v_eci(:,1)] ;

[x0_Tango,resnorm_Tango,residual_Tango,exitflag_Tango,~,~,jacobian_Tango] = ...
lsqnonlin(@costFcn, x0_guess, [],[],options,data,mu,fun,Mango.initial_epoch,ode_opt,J2,Re);

sol4 = x0_guess - x0_Tango;
P_ls_Tango = resnorm_Tango/(length(residual_Tango)-length(x0_Tango))...
    .*inv(jacobian_Tango.'*jacobian_Tango);
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

x = xx(1);
y = xx(2);
z = xx(3);
r_norm = norm([x;y;z]);
% RHS
dxx(1:3) = xx(4:6);
dxx(4:6) = -mu/r_norm^3 * [x;y;z];

dxx = dxx';

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

function [sat_coord,varargout] = scLocalCoordinates(station,states)

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
if nargout == 2
    measures = mvnrnd(sat_coord(1:3,:)',station.R);
    measures = measures(measures(:,3)>=station.min_el,:);
    vis_time = station.visibility_time(measures(:,3)>=station.min_el)';
    varargout{1} = [vis_time, measures];
end

end

function [visibility_time, coord] = visibility(sat_coord,min_el,t_span)

    visibility_condition = sat_coord(3,:)>=min_el; 
    visibility_time = t_span(visibility_condition);
    coord = sat_coord(:,visibility_condition);
end

function [station1,station2,sat,station1_measure,station2_measure] = ...
    simulateMeasurements(TLE,station1,station2)

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
sat.time = [station1.visibility_time, station2.visibility_time];
sat.time_TLE = (sat.time - sat.initial_epoch)./60;

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
ttt = cspice_unitim(sat.time(i), 'ET', 'TDT')/cspice_jyear()/100;

% rotation from TEME to ECI
[sat.r_eci(:,i),sat.v_eci(:,i),~] = teme2eci(sat.r_teme(:,i),sat.v_teme(:,i), ateme, ttt, dpsi,deps);
end

% simulate measurements

station1.visibility_length = length(station1.visibility_time);
station2.visibility_length = length(station2.visibility_time);

[station1.eci,station1.eci2topo] = stationCoordinates(station1,station1.visibility_time);
coord = [sat.r_eci(:,1:station1.visibility_length);sat.v_eci(:,1:station1.visibility_length)];
[~,station1_measure] =...
    scLocalCoordinates(station1,coord);

[station2.eci,station2.eci2topo] = stationCoordinates(station2,station2.visibility_time);
coord = [sat.r_eci(:,end-station2.visibility_length+1:end);sat.v_eci(:,end-station2.visibility_length+1:end)];
[~,station2_measure] =...
    scLocalCoordinates(station2,coord);
end

function residual = costFcn (x0,stations,mu,fun,t0,ode_opt,varargin)

if nargin == 8
    J2 = varargin{1};
    Re = varargin{2};
end
t_span = zeros(10000,1);
measure = zeros(10000,3);
weight = zeros(3,3,10000);
station_eci = zeros(6,10000);
station_eci2topo = zeros(6,6,10000);
index_old = 0;

for i = 1:size(stations,1)
    index_new = index_old + size(stations{i,2},1);
    t_span(index_old+1:index_new) = stations{i,2}(:,1);
    measure(index_old+1:index_new,:) = stations{i,2}(:,2:4);
    weight(:,:,index_old+1:index_new) = stations{i,1}.W_m.*ones(3,3,index_new-index_old);
   
    [station_eci(:,index_old+1:index_new),station_eci2topo(:,:,index_old+1:index_new)] =...
        stationCoordinates (stations{i,1}, stations{i,2}(:,1)');

    index_old = index_new;
end
% get sorted measurements

t_span = t_span(1:index_new);
[t_span,index] = sort(t_span);

measure = measure(index,:);
weight = weight(:,:,index);
station.eci = station_eci(:,index);
station.eci2topo = station_eci2topo(:,:,index);

% propagate orbit
t_ode = [t0;t_span];
[~,x] = ode113(fun,t_ode,x0,ode_opt);
x = x(2:end,:);
% convert model to range,Az,El
model_coord = scLocalCoordinates(station,x');
% compute residual
residual = zeros(length(t_span),3);
for i = 1:length(t_span)
 residual(i,:) = weight(:,:,i)*[model_coord(1,i)'-measure(i,1),...
     angdiff(model_coord(2:3,i)'*pi/180,measure(i,2:3)*pi/180) ]';
end

end