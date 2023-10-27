% Spacecraft Guidance and Navigation (2023/24)
% Assignment 1, exercise 2
% Author: Giuseppe Brentino

clearvars; close all; clc;

addpath('.\kernels')
addpath('.\mice\src\mice')
addpath('.\mice\lib')

cspice_furnsh('kernels\naif0012.tls'); % (LSK)
cspice_furnsh('kernels\de432s.bsp');   % (SPK)
cspice_furnsh('kernels\20099942_Apophis.bsp'); % (PCK)
cspice_furnsh('kernels\gm_de432.tpc');  
cspice_furnsh('kernels\pck00010.tpc'); 
%% Ex 1

% define close approach time window
close_approach_epoch = datetime(2029,1,1):hours(1):datetime(2029,7,31);

% compute Apophis distance with Moon-Eart-Sun during the defined time
% window

% a)
    date =  cspice_str2et( char(close_approach_epoch ) );
    earth = cspice_spkpos('20099942',date,'IAU_EARTH','NONE','EARTH');
    moon = cspice_spkpos('20099942',date,'IAU_EARTH','NONE','MOON');
    sun = cspice_spkpos('20099942',date,'IAU_EARTH','NONE','SUN');
    
    apophis_earth = vecnorm(earth,2,1);
    apophis_moon = vecnorm(moon,2,1);
    apophis_sun = vecnorm(sun,2,1);

% b)
    earth_apophis_sun_angle = acosd(dot(earth,sun)./(apophis_sun.*apophis_earth ));

% c)
min_time = fmincon(@apophisEarthFindMin,date(round(length(date)/2)),[],[],[],[],date(1),date(end),[]);
time_window = (min_time-6*3600):(min_time+6*3600);
rec_coord = cspice_spkpos('20099942',time_window,'IAU_EARTH','NONE','EARTH');
%1) Call cspice_bodvrd to retrieve Earth's radii:
Re=cspice_bodvrd('399','RADII',3);
%2) Compute flatness f=(Re-Rp)/Re:
fe=(Re(1)-Re(3))/Re(1);

% doesn't take into account earth rotation, DA VEDERE
[lon,lat,~] =  cspice_recgeo( rec_coord, Re(1), fe );


% plots
figure()

subplot(1,3,1)
hold on
grid on
plot(date,apophis_sun)
title('Sun-Apophis distance')
xlabel('Time [s]')
ylabel('Distance [km]')

subplot(1,3,2)
hold on
grid on
plot(date,apophis_moon)
title('Moon-Apophis distance')
xlabel('Time [s]')
ylabel('Distance [km]')

subplot(1,3,3)
hold on
grid on
plot(date,apophis_earth)
title('Earth-Apophis distance')
xlabel('Time [s]')
ylabel('Distance [km]')

figure()
hold on
grid on
plot(date,earth_apophis_sun_angle)
title('earth-apophis-sun-angle')

figure()
hold on
grid on
plot(rad2deg(lon),rad2deg(lat))
%% Ex 2 
clearvars; close all; clc;

close_approach_epoch = datetime(2029,1,1):caldays(1):datetime(2030,1,1);
date =  cspice_str2et( char(close_approach_epoch ) );

x0 = cspice_spkezr('20099942',date(1),'J2000','NONE','SUN');

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
% x = state1, state2, state3, t0,tdsm,timp;
LWO = cspice_str2et( char(datetime(2024,10,1)) );
LWC = cspice_str2et( char(datetime(2025,2,1)) );

DSMO = cspice_str2et( char(datetime(2024,10,1)+calmonths(6)) );
DSMC = cspice_str2et( char(datetime(2024,10,1)+calmonths(18)) );

IMPO = cspice_str2et( char(datetime(2028,8,1)) );
IMPC = cspice_str2et( char(datetime(2029,2,28)) );

LW_vect = linspace(LWO,LWC,300);
earth_states = cspice_spkezr('EARTH',LW_vect,frame,'NONE',center);
min_earth_states = min(earth_states,[],2);
max_earth_states = max(earth_states,[],2);

IMP_vect = linspace(IMPO,IMPC,300);
apophis_states = cspice_spkezr('20099942',LW_vect,frame,'NONE',center);
min_apo_states = min(apophis_states,[],2);
max_apo_states = max(apophis_states,[],2);

x1_lb = min([min_earth_states min_apo_states],[],2) - [1e3*ones(3,1);10*ones(3,1)];
x1_ub = max([max_earth_states max_apo_states],[],2) + [1e3*ones(3,1);10*ones(3,1)];

lb = [x1_lb; x1_lb; x1_lb; LWO; DSMO; IMPO];
ub = [x1_ub; x1_ub; x1_ub; LWC; DSMC; IMPC];

%Initial conditions
x0 = zeros(21,1);
x0(19) = (LWO+LWC)/2;
x0(20) = (DSMO+DSMC)/2;
x0(21) = (IMPO+IMPC)/2;
x0(1:6) = cspice_spkezr('Earth', x0(19), frame, 'NONE', center);
ode_opt = odeset('RelTol',1e-11,'AbsTol',1e-12);
[~,x2_guess] = ode78(@scPropagator,[x0(19) x0(20)],x0(1:6),ode_opt,mu);
x0(7:12) = x2_guess(end,:)';
[~,x3_guess] = ode78(@scPropagator,[x0(20) x0(21)],x0(7:12),ode_opt,mu);
x0(13:18) = x3_guess(end,:)'; 

opt = optimoptions('fmincon','Display','iter-detailed','UseParallel',true);
parfor i = 1:4
    cspice_furnsh('kernels\naif0012.tls'); % (LSK)
    cspice_furnsh('kernels\de432s.bsp');   % (SPK)
    cspice_furnsh('kernels\20099942_Apophis.bsp'); % (PCK)
    cspice_furnsh('kernels\gm_de432.tpc');
    cspice_furnsh('kernels\pck00010.tpc');
end
[states,fval] = fmincon(@(x)guidance(x,mu,frame,center,bodies),x0,[],[],[],[],lb,ub,@(x)nonlcon(x,mu),opt);

%% functions Ex 1

function [distance] = apophisEarthFindMin(t)    
    distance = cspice_spkpos('20099942',t,'IAU_EARTH','NONE','EARTH');
    distance = vecnorm(distance);
end

%% functions Ex 2

function dx = scPropagator(t,x,mu)
    
    r = x(1:3);
    v = x(4:6);
    r_norm = norm(r);
    
    dx(1:3) = v;
    dx(4:6) = - mu/r_norm^3 * r;

    dx = dx';

end

function [dx] = nbody_rhs(t, x, bodies, frame)
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

bodies = cell(size(labels));

for i = 1:length(labels)
    bodies{i}.name = labels{i};
    bodies{i}.GM = cspice_bodvrd(labels{i}, 'GM',1);
end

end

function J = guidance(x,mu,frame,center,bodies)

t_dsm = x(20);
t_imp = x(21);

x2 = x(7:12);
ode_opt = odeset('RelTol',1e-11,'AbsTol',1e-12);
[~,x_sc_2] = ode78(@scPropagator,[t_dsm t_imp],x2,ode_opt,mu);

x_Ai = cspice_spkezr('20099942',t_imp,frame,'NONE',center) + [zeros(3,1); x_sc_2(end,4:6)'*5e-5];
opt = odeset('RelTol',1e-8,'AbsTol',1e-9,'Events',@(t,s) minDistanceEvent(t,s));
[T,s] = ode78(@(t,s) nbody_rhs(t,s,bodies,frame),[t_imp t_imp*10],x_Ai,opt);

J = -norm( s(end,1:3) - cspice_spkpos('EARTH',T(end),frame,'NONE',center)' );  

end

function [C, Ceq] = nonlcon(x,mu)
    
    s_e = cspice_spkezr('EARTH',x(19),'ECLIPJ2000','NONE','SSB');
    Ceq(1:3) = x(1:3) - s_e(1:3) ;
   
    ode_opt = odeset('RelTol',1e-11,'AbsTol',1e-12);
    [~,x_sc_1] = ode78(@scPropagator,[x(19) x(20)],x(1:6),ode_opt,mu);
    Ceq(4:6) = x(7:9) -x_sc_1(end,1:3)';

    ode_opt = odeset('RelTol',1e-11,'AbsTol',1e-12);
    [~,x_sc_2] = ode78(@scPropagator,[x(20) x(21)],x(7:12),ode_opt,mu);
    Ceq(7:12) = x(13:18) - x_sc_2(end,:)';

    Ceq(13:15) = x(13:15) - cspice_spkpos('20099942',x(21),'ECLIPJ2000','NONE','SSB');

    C = norm(x(4:6)-s_e(4:6)) + norm(x_sc_1(end,4:6)'-x(10:12)) - 5;

end

function [value,isterminal,direction] = minDistanceEvent(t,x)
x_e = cspice_spkezr('EARTH',t,'ECLIPJ2000','NONE','SSB');
value = dot(x_e(1:3)-x(1:3),x_e(4:6)-x(4:6));  
isterminal = 1;         % stop at local minimum
direction  = 1;         % [local minimum, local maximum]
end
