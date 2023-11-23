% Spacecraft Guidance and Navigation (2023/24)
% Assignment 1, exercise 2
% Author: Giuseppe Brentino

clearvars; close all; clc;
rng default;

% load kernels
addpath('.\kernels')
addpath('.\mice\src\mice')
addpath('.\mice\lib')

cspice_furnsh('kernels\naif0012.tls'); % (LSK)
cspice_furnsh('kernels\de432s.bsp');   % (SPK)
cspice_furnsh('kernels\20099942_Apophis.bsp'); % (PCK)
cspice_furnsh('kernels\gm_de432.tpc');
cspice_furnsh('kernels\pck00010.tpc');

% set plot properties
plotStyle;

%% Ex 1
clearvars; close all; clc;

% define close approach time window
close_approach_epoch = datetime(2029,1,1):hours(1):datetime(2029,7,31);

% compute earth,moon and sun relative position wrt apophis
date =  cspice_str2et( char(close_approach_epoch ) );
earth = cspice_spkpos('20099942',date,'IAU_EARTH','NONE','EARTH');
moon = cspice_spkpos('20099942',date,'IAU_EARTH','NONE','MOON');
sun = cspice_spkpos('20099942',date,'IAU_EARTH','NONE','SUN');
 
% compute Apophis distance with Moon-Eart-Sun during the defined time
% window
apophis_earth = vecnorm(earth,2,1);
apophis_moon = vecnorm(moon,2,1);
apophis_sun = vecnorm(sun,2,1);

% compute earth-apophis-sun angle
earth_apophis_sun_angle = acosd(dot(earth,sun)./(apophis_sun.*apophis_earth ));

% find minimum distance epoch
min_time = fmincon(@apophisEarthFindMin,date(round(length(date)/2)),[],[],[],[],date(1),date(end),[]);

% compute groundtrack
time_window = (min_time-6*3600):(min_time+6*3600);
rec_coord = cspice_spkpos('20099942',time_window,'IAU_EARTH','NONE','EARTH');
Re=cspice_bodvrd('399','RADII',3);
fe=(Re(1)-Re(3))/Re(1);
[lon,lat,~] = cspice_recgeo( rec_coord, Re(1), fe );

%%% Plots
plot_date_vec = linspace(date(1),date(end),3);
plot_dates=cspice_et2utc(plot_date_vec,'C', 1e-3);

% plot distances
figure()

subplot(1,3,1)
hold on
grid on
plot(date,apophis_sun)
title('Sun-Apophis distance','HorizontalAlignment','left')
xlabel('Date')
ylabel('Distance [km]')
xticks(plot_date_vec)
xticklabels(plot_dates(:,1:11))

subplot(1,3,2)
hold on
grid on
plot(date,apophis_moon)
title('Moon-Apophis distance','HorizontalAlignment','left')
xlabel('Date')
ylabel('Distance [km]')
xticks(plot_date_vec)
xticklabels(plot_dates(:,1:11))

subplot(1,3,3)
hold on
grid on
plot(date,apophis_earth)
title('Earth-Apophis distance','HorizontalAlignment','left')
xlabel('Date')
ylabel('Distance [km]')
xticks(plot_date_vec)
xticklabels(plot_dates(:,1:11))

% plot angle
figure()
hold on
grid on
plot(date,earth_apophis_sun_angle)
xlabel('Date')
ylabel('Angle [deg]')
xticks(plot_date_vec)
xticklabels(plot_dates(:,1:11))

% plot groundtrack
lon = rad2deg(lon);
lat = rad2deg(lat);

figure()
hold on
photo = imread('Earth.png');
image('XData', [-180, 180], 'YData', [-90, 90], 'CData', flipud(photo));

% remove horizontal line in the plot when longitude changes sign
for j = 2:length(lon)
    if lon(j-1)*lon(j) < 0
        lon = [lon(1:j-1) NaN lon(j:end)];
        lat = [lat(1:j-1) NaN lat(j:end)];
    end
end

% compute coordinates at closest approach
CA_coord = cspice_spkpos('20099942',min_time,'IAU_EARTH','NONE','EARTH');
[lon_CA,lat_CA,~] =  cspice_recgeo( CA_coord, Re(1), fe );

plot(lon,lat,'r','DisplayName','Groundtrack')
plot(rad2deg(lon_CA),rad2deg(lat_CA),'o','MarkerFaceColor','#EDB120',...
    'MarkerEdgeColor','#EDB120','DisplayName','Closest approach')
legend('Location','northoutside');

ylabel('Latitude [deg]')
xlabel('Longitude [deg]')
axis equal
xlim([-180,180])
ylim([-90,90])
%% Ex 2
clearvars -EXCEPT min_time; close all; clc;

% compute closest approach epoch
close_approach_epoch = datetime(2029,1,1):caldays(1):datetime(2030,1,1);
date =  cspice_str2et( char(close_approach_epoch ) );


labels = {'Sun';
    'Mercury';
    'Venus';
    'Earth';
    'Moon';
    'Mars Barycenter';
    'Jupiter Barycenter';
    'Saturn Barycenter';
    'Uranus Barycenter';
    'Neptune Barycenter';
    'Pluto Barycenter'};

% Initialize propagation data (same as regular n-body)
bodies = nbody_init(labels);
mu = cspice_bodvrd('SUN', 'GM',1);

% select integration frame string (SPICE naming convention)
frame = 'ECLIPJ2000';
center = 'SSB';

% Time windows upper and lower bounds
LWO = cspice_str2et( char(datetime(2024,10,1)) );
LWC = cspice_str2et( char(datetime(2025,2,1)) );

DSMO = cspice_str2et( char(datetime(2024,10,1)+calmonths(6)) );
DSMC = cspice_str2et( char(datetime(2025,2,1)+calmonths(18)) );

IMPO = cspice_str2et( char(datetime(2028,8,1)) );
IMPC = cspice_str2et( char(datetime(2029,2,28)) );

% max and minimum earth states
LW_vect = linspace(LWO,LWC,300);
earth_states = cspice_spkezr('EARTH',LW_vect,frame,'NONE',center);
min_earth_states = min(earth_states,[],2);
max_earth_states = max(earth_states,[],2);

% max and minimum Apophis states
IMP_vect = linspace(IMPO,IMPC,300);
apophis_states = cspice_spkezr('20099942',LW_vect,frame,'NONE',center);
min_apo_states = min(apophis_states,[],2);
max_apo_states = max(apophis_states,[],2);

% compute position and velocity bounds as the minimum and maximum of the
% departure and target body plus an offset
x1_lb = min([min_earth_states min_apo_states],[],2) - [1e3*ones(3,1);10*ones(3,1)];
x1_ub = max([max_earth_states max_apo_states],[],2) + [1e3*ones(3,1);10*ones(3,1)];

lb = [x1_lb; x1_lb; x1_lb; LWO; DSMO; IMPO];
ub = [x1_ub; x1_ub; x1_ub; LWC; DSMC; IMPC];

% load kernels in each core to allow fmincon to compute gradients in
% parallel
num_workers = parcluster('processes').NumWorkers;
parfor i = 1:num_workers
    cspice_furnsh('kernels\naif0012.tls'); % (LSK)
    cspice_furnsh('kernels\de432s.bsp');   % (SPK)
    cspice_furnsh('kernels\20099942_Apophis.bsp'); % (PCK)
    cspice_furnsh('kernels\gm_de432.tpc');
    cspice_furnsh('kernels\pck00010.tpc');
end

% fmincon and ODE options
opt = optimoptions('fmincon','Display','iter-detailed','UseParallel',true);
ode_opt = odeset('RelTol',1e-11,'AbsTol',1e-12);

% initialize first guess
x0 = zeros(21,1);

% compute initial guesses
x0(19) =   (LWO+LWC)/2 - (LWC-LWO)*0.4;
x0(20) =   (DSMO+DSMC)/2 - (LWC-LWO)*0.4;
x0(21) =  (IMPO+IMPC)/2 - (LWC-LWO)*0.4;
x0(1:6) = cspice_spkezr('Earth', x0(19), frame, 'NONE', center);

[~,x2_guess] = ode78(@scPropagator,[x0(19) x0(20)],x0(1:6),ode_opt,mu);
x0(7:12) = x2_guess(end,:)';
[~,x3_guess] = ode78(@scPropagator,[x0(20) x0(21)],x0(7:12),ode_opt,mu);
x0(13:18) = x3_guess(end,:)';

% optimization

[states,fval,exitflag] = fmincon(@(x)guidance(x,frame,center,bodies),x0,[],[],[],[],lb,ub,@(x)nonlcon(x,mu),opt);

%%% retrive data

% first delta v manoeuvre
x_e =  cspice_spkezr('Earth',states(19),frame,'NONE',center);
dv0 = states(4:6) - x_e(4:6);
dv0_norm = norm(dv0);

% second delta v manoeuvre
[~,x_sc_prop] = ode78(@scPropagator,[states(19) states(20)],states(1:6),ode_opt,mu);
dv_dsm = states(10:12)-x_sc_prop(end,4:6)';
dv_dsm_norm = norm(dv_dsm);

% relevant times
x_Ai = cspice_spkezr('20099942',states(21),frame,'NONE',center) + [zeros(3,1); states(16:18)*5e-5];
opt = odeset('RelTol',1e-8,'AbsTol',1e-9,'Events',@(t,s) minDistanceEvent(t,s));
[T_Aimp,A_aimp] = ode78(@(t,s) nbody_rhs(t,s,bodies,frame),[states(21) states(21)*10],x_Ai,opt);
tdca = cspice_et2utc(T_Aimp(end),'C',1);
t0 = cspice_et2utc(states(19),'C',1);
tdsm = cspice_et2utc(states(20),'C',1);
timp = cspice_et2utc(states(21),'C',1);

% earth and apophis orbits during the mission
time_vec = linspace(states(19),min_time,500);
time_vec_a = linspace(states(19),states(21),500);
earth = cspice_spkpos('EARTH',time_vec,frame,'NONE',center);
apo = cspice_spkpos('20099942',time_vec_a,frame,'NONE',center);
apo = [apo A_aimp(:,1:3)'];

% plot 3d orbits
sc_i = x_sc_prop(:,1:3)';
[~,x_sc_prop2] = ode78(@scPropagator,[states(20) states(21)],states(7:12),ode_opt,mu);
sc = [sc_i x_sc_prop2(:,1:3)'];

figure()
hold on
grid on
plot3(apo(1,:),apo(2,:),apo(3,:),'Color',[0.2,0.2,0.2])
plot3(earth(1,:),earth(2,:),earth(3,:),'Color','b')
plot3(sc(1,:),sc(2,:),sc(3,:),'Color','r')
plot3(sc(1,1),sc(2,1),sc(3,1),'o','MarkerSize',5,'MarkerFaceColor','g')
plot3(sc_i(1,end),sc_i(2,end),sc_i(3,end),'o','MarkerSize',5,'MarkerFaceColor','c')
plot3(sc(1,end),sc(2,end),sc(3,end),'o','MarkerSize',5,'MarkerFaceColor','m')
axis equal

% plot distance
figure()
hold on
grid on
plot_date_vec = linspace(time_vec(1),T_Aimp(end),5);
plot_dates=cspice_et2utc(plot_date_vec,'C', 1e-3);

Re=cspice_bodvrd('399','RADII',3);
dreal = vecnorm(cspice_spkpos('20099942',time_vec,'IAU_EARTH','NONE','EARTH'),2,1)./Re(1);
dmod = vecnorm(cspice_spkpos('EARTH',T_Aimp',frame,'NONE',center)-A_aimp(:,1:3)',2,1)./Re(1);
plot(time_vec,dreal,'b-','LineWidth',2)
plot(T_Aimp,dmod,'r--','LineWidth',2)
xticks(plot_date_vec)
xticklabels(plot_dates(:,1:11))
ylabel('Distance [Earth Radii]')
DCA = min(dmod);

axes('position',[0.61 0.61 0.3 0.3])
plot_date_vec = linspace(time_vec(1),T_Aimp(end),7000);
plot_dates=cspice_et2utc(plot_date_vec,'C', 1e-3);

box on % put box around new pair of axes
hold on
grid on
plot(time_vec,dreal,'b-','LineWidth',2)
plot(T_Aimp,dmod,'r--','LineWidth',2)
ylim([5 50])
xlim([T_Aimp(find(dmod<50,1,'first')) T_Aimp(end)])
xticks(plot_date_vec)
xticklabels(plot_dates(:,6:20))
set(gca,'color','w');
ylabel('Distance [Earth Radii]')

%% functions Ex 1

function [distance] = apophisEarthFindMin(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compute the distance between the Earth and Apophis at time
% t.
%
% INPUT: 
% t[1]: scalar, epoch at which the function computes the distance
%
% OUTPUT: 
% distance[1]: distance between the Earth and Apophis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance = cspice_spkpos('20099942',t,'IAU_EARTH','NONE','EARTH');
distance = norm(distance);

end

%% functions Ex 2

function dx = scPropagator(~,x,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the derivative of the state vector for the 
% two-body problem, which is used as the input to an ODE solver.
%
% INPUTS:
% - ~: Placeholder for the time variable (not used in this function).
% - x[6x1]: State vector containing position and velocity components.
%             [r1; r2; r3; v1; v2; v3]
% - mu[1]: Gravitational parameter of the central body.
%
% OUTPUT:
% - dx[6x1]: Derivative of the state vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define states
r = x(1:3);
v = x(4:6);
r_norm = norm(r);

% compute derivatives
dx(1:3) = v;
dx(4:6) = - mu/r_norm^3 * r;

dx = dx';

end

function [dx] = nbody_rhs(t, x, bodies, frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the right-hand side of the ordinary differential 
% equation (ODE) for the n-body problem.
%
% INPUTS:
% - t[1]: Time variable.
% - x[6x1]: State vector containing position and velocity components.
% - bodies[cell array]: Cell array containing information about celestial 
%                      bodies (name, gravitational parameter GM).
% - frame[char]: Reference frame for computing positions (e.g., 'ECLIPJ2000').
%
% OUTPUT:
% - dx[6x1]: Derivative of the state vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize right-hand-side
dx = zeros(6,1);
% Extract the object position from state x
rr_ssb_obj = x(1:3);
%r = norm(rr);
% Loop over all bodies
for i = 1:length(bodies)

    % Retrieve position and velocity of i-th celestial body wrt Solar
    % System Barycentre in inertial frame
    %bodies{i}.name
    rv_ssb_body = cspice_spkezr(bodies{i}.name,t,frame,'NONE','SSB');

    % Extract object position wrt. i-th celestial body
    rr_body_obj = rr_ssb_obj - rv_ssb_body(1:3);
    % Compute square distance and distance
    %r_i = norm(rr_i);
    dist2 = dot(rr_body_obj, rr_body_obj);
    dist = sqrt(dist2);
    % Compute the gravitational acceleration using Newton's law
    % GM = bodies{i}.GM;
    aa_grav = - bodies{i}.GM*rr_body_obj / (dist*dist2);

    % Position derivative is object's velocity
    dx(1:3) = x(4:6);
    % Sum up acceleration to right-hand-side
    dx(4:6) = dx(4:6) + aa_grav;

end
end

function [bodies] = nbody_init(labels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes celestial bodies for the n-body problem.
%
% INPUT:
% - labels[cell array]: Cell array containing names of celestial bodies.
%
% OUTPUT:
% - bodies[cell array]: Cell array containing information about celestial 
%                      bodies (name, gravitational parameter GM).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bodies = cell(size(labels));

for i = 1:length(labels)
    bodies{i}.name = labels{i};
    bodies{i}.GM = cspice_bodvrd(labels{i}, 'GM',1);
end

end

function J = guidance(x,frame,center,bodies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the guidance objective function for trajectory optimization.
%
% INPUTS:
% - x[21x1]: Vector containing optimization parameters.
% - frame[char]: Reference frame for computing positions (e.g., 'ECLIPJ2000').
% - center[char]: Central body for position computation (e.g., 'SSB').
% - bodies[cell array]: Cell array containing information about celestial 
%                      bodies (name, gravitational parameter GM).
%
% OUTPUT:
% - J[1]: Guidance objective function value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time of impact
t_imp = x(21);

% compute the new Apophis state vector after impact
x_Ai = cspice_spkezr('20099942',t_imp,frame,'NONE',center) + [zeros(3,1); x(16:18)*5e-5];

% Define options for the ODE solver
opt = odeset('RelTol',1e-8,'AbsTol',1e-9,'Events',@(t,s) minDistanceEvent(t,s));
[T,s] = ode78(@(t,s) nbody_rhs(t,s,bodies,frame),[t_imp t_imp*10],x_Ai,opt);

% Compute the guidance objective function as the negative of the final
% distance between Apophis and the Earth
J = -norm( s(end,1:3) - cspice_spkpos('EARTH',T(end),frame,'NONE',center)' );

end

function [C, Ceq] = nonlcon(x,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear constraint function for trajectory optimization.
%
% INPUTS:
% - x[21x1]: Vector containing optimization parameters.
% - mu[1]: Gravitational parameter of the central body.
%
% OUTPUTS:
% - C[1]: Inequality constraint value.
% - Ceq[15x1]: Equality constraint values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize parameters
t0 = x(19);
tdsm = x(20);
timp = x(21);
x1 = x(1:6);
x2 = x(7:12);
x3 = x(13:18);
Ceq = zeros(15,1);

% Compute position of the Earth at t0
s_e = cspice_spkezr('EARTH',t0,'ECLIPJ2000','NONE','SSB');

% Enforce that the initial position of the asteroid (x(1:3)) is equal to
% the position of the Earth (s_e(1:3)) at the initial time t0
Ceq(1:3) = x(1:3) - s_e(1:3) ; 

% Propagate spacecraft dynamics from t0 to tdsm
ode_opt = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,x_sc_1] = ode78(@scPropagator,[t0 tdsm],x1,ode_opt,mu);

% Enforce that the final position of the spacecraft  (x_sc_1(end, 1:3))
% after propagating from t0 to tdsm is equal to the initial position of the
% second segment (x2(1:3))
Ceq(4:6) = x_sc_1(end,1:3)' - x2(1:3);

% Propagate spacecraft dynamics from tdsm to timp
[~,x_sc_2] = ode78(@scPropagator,[tdsm timp],x2,ode_opt,mu);

% Enforce that the final state of the spacecraft (x_sc_2(end, :)) after
% propagating from tdsm to timp is equal to the initial state of the third
% segment (x3)
Ceq(7:12) = x_sc_2(end,:)'-x3;

% Enforce that the final position of the third segment (x3(1:3)) is equal 
% to the position of Apophis at time timp
Ceq(13:15) = x3(1:3) - cspice_spkpos('20099942',timp,'ECLIPJ2000','NONE','SSB');

% Enforce a maximum total delta V of 5 km/s
C = norm(x(4:6)-s_e(4:6)) + norm(x_sc_1(end,4:6)'-x(10:12)) - 5;

end

function [value,isterminal,direction] = minDistanceEvent(t,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines an event to monitor the minimum distance between
% the Earth and another body orbiting the solar system barycenter.
%
% INPUTS:
% - t[1]: Time variable.
% - x[6x1]: State vector of the body.
%
% OUTPUTS:
% - value[1]: The dot product of the vectors connecting the spacecraft to 
%             the Earth and the relative velocity: it's the first
%             derivative of the relative position.
% - isterminal[1]: Indicates whether the integration should terminate when
%                  the event is triggered (1 for true, 0 for false).
% - direction[1]: Specifies whether the event should be detected in the
%                 increasing (1) or decreasing (-1) direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_e = cspice_spkezr('EARTH',t,'ECLIPJ2000','NONE','SSB');
value = dot(x_e(1:3)-x(1:3),x_e(4:6)-x(4:6));
isterminal = 1;         % stop at local minimum
direction  = 1;         % [local minimum, local maximum]
end

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










