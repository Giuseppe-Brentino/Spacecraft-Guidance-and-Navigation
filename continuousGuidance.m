% Spacecraft Guidance and Navigation (2023/24)
% Assignment 1, exercise 3
% Author: Giuseppe Brentino

%% Ex 2 
clearvars; close all; clc;
plotStyle;

addpath('.\kernels')
addpath('.\mice\src\mice')
addpath('.\mice\lib')

cspice_furnsh('kernels\naif0012.tls'); % (LSK)
cspice_furnsh('kernels\de432s.bsp');   % (SPK)
cspice_furnsh('kernels\gm_de432.tpc');  

% set Launch epoch
t_launch = cspice_str2et('May 28 14:13:09.000 UTC 2023');

% define data
m0 = 1000;          % Spacecraft mass [kg]
Tmax = 800e-6;      % Maximum thrust available [kN]   
Isp = 3120;         % Specific Impulse [s]
g0 = 9.81e-3;       % Standard acceleration of gravity [km/s^2]
data.center = 'SUN';
data.frame = 'ECLIPJ2000'; 

% Adimentionalization units
mu_s = cspice_bodvrd('SUN','GM',1);
l = cspice_convrt(1,'AU','km');
m = m0;
t = sqrt(l^3/(mu_s));

% set initial conditions of the spacecraft
xi = cspice_spkezr('Earth',t_launch,data.frame,'NONE',data.center);

x0 = [xi;m0];
x0 = [x0(1:3)/l;x0(4:6)*(t/l);x0(7)/m];

% Adimentionalize data
data.Tmax = Tmax*((t^2)/(m*l));
data.Isp = Isp/t;
data.m0 = m0/m;
data.g0 = g0*(t^2)/l;
data.mu_s = mu_s*(t^2)/(l^3);
data.t0 = t_launch/t;
data.l = l;
data.m = m;
data.t = t;
data.ode_opt = odeset('RelTol',1e-12,'AbsTol',1e-13);

%% Ex 3
close all; clearvars -except data x0; clc;
rng default
N_iter = 10;

% find correct costate initial conditions and propagate the s/c orbit
[sol,errors,T,xx,venus,Timedata] = computeOrbit(N_iter,data,x0);
   
 % plot s/c trajectory
plotOrbits(venus,data,xx);

 % check hamiltonian
 H = zeros(length(T),1);
 for i = 1:length(T)
     xH0 = xx(i,:);
     dyn = TPBVP(T(i),xH0,data);
     H(i) = 1 + dot(xx(i,8:end),dyn(1:7));
 end
 figure()
 hold on
 grid on
 plot(T,(H-H(1))./H(1))
 xlabel('Time of flight')
 ylabel('Relative error')

%% Ex 4 
% close all; clc; clearvars -except data x0 sol;
rng default;

% Maximum thrust vector, from 700 to the required 500mN to exploit
% numerical continuation
Tmax_vect = linspace(700,500,5)*1e-6;

% Find initial costates for every thrust level of Tmax to exploit numerical
% continuation
N_iter = 10;

for j = 1:length(Tmax_vect)
   
    %adimentionalize thrust
    data.Tmax_vect = Tmax_vect(j)*((data.t^2)/(data.m*data.l));

    % new initial guess is the solution of the previous iteration
    guess = sol;

    % find correct costate initial conditions and propagate the s/c orbit
    [sol,errors,T,xx,venus,Timedata] = computeOrbit(N_iter,data,x0,guess);

end

% plot trajectory
plotOrbits(venus,data,xx);

%% functions

function dx = TPBVP(~,x,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function integrates the Two Point Boundary Value Problem for the
% specific case of a time-optimal, continuous-thrust spacecraft with the
% objective of matching the position and velocity of another celestial body
% at the final time. The Two-body problem simplification of the dynamics is
% employed.
%
% INPUT:
% - t[1]: time instant. It must be the first input in order to use
%         this function inside an ode solver, but its not required in this
%         case.
% - x[14,1]: State vector.
% - data[struct]: Struct containing adimentionalization constants and other
%                 spacecraft and environment parameters.
%
% OUTPUT:
% -dx[42,1]: time derivative of the state vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract data
mu = data.mu_s;
Isp = data.Isp;
g0 = data.g0;
Tmax = data.Tmax;

rr = x(1:3);
vv = x(4:6);
m = x(7);
ll_r = x(8:10);
ll_v = x(11:13);

r = norm(rr);
l_v = norm(ll_v);

% compute RHS
dx = zeros(14,1);

dx(1:3) = vv;
dx(4:6) = -mu/r^3*rr -Tmax/m * ll_v/l_v;
dx(7) = -Tmax/(Isp*g0);
dx(8:10) = -3*mu/r^5 * dot(rr,ll_v)*rr + mu/r^3*ll_v;
dx(11:13) = -ll_r;
dx(14) = -l_v*Tmax/m^2;

end

function Fun = optimization(guess,x0,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used inside a zero-finding solver in order to find the
% initial costates that satisfy the requirements of this specific problem,
% which are:
% - spacecraft final position and velocity must match the ones of
%   venus at the final time;
% - the mass co-state must be equal to zero at the final time;
% - the transversality condition must be satisfied.
% 
% INPUT:
% - guess[8,1]: Initial guess of the costates and final time.
% - x0[7,1]: Initial state vector of the spacecraft.
% - data[struct]: Struct containing adimentionalization constants and other
%                 spacecraft and environment parameters.
%
% OUTPUT: 
% Fun[8,1]: vector containing the functions that must be equal to 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract data
mu = data.mu_s;
Isp = data.Isp;
g0 = data.g0;
Tmax = data.Tmax;
t0 = data.t0;
frame = data.frame;
center = data.center;
TU = data.t;
LU = data.l;

% propagate s/c dynamics
s0 = [x0;guess(1:7)];
tf = guess(8);
[~,states] = ode113(@TPBVP,[t0 tf],s0,data.ode_opt,data);

% extract final conditions
venus = cspice_spkezr('VENUS',tf*TU,frame,'NONE',center);
venus = venus./[LU*ones(3,1);(LU/TU)*ones(3,1)]; 
rr_f = states(end,1:3)';
vv_f = states(end,4:6)';
ll_f = states(end,8:14)';
r_f = norm(rr_f);
lv_f = norm(ll_f(4:6));
m_f = states(end,7);

%Hamiltonian
   dyn       =  [vv_f;
                -mu/r_f^3*rr_f-Tmax/m_f * ll_f(4:6)/lv_f;
                -Tmax/(Isp*g0)];
   H = 1 + dot(ll_f,dyn);

% final venus conditions 
psi_dot = [venus(4:6);-mu/norm( venus(1:3) )^3 * venus(1:3)];

% final function to evaluate
Fun = [ rr_f-venus(1:3);
        vv_f-venus(4:6);
        ll_f(7);
        H - dot(ll_f(1:6),psi_dot);
    ];

end

function [sol,errors,T,xx,venus,Timedata] = computeOrbit(N_iter,data,x0,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the correct initial costates and final time that
% satisfy the requirements, also returning the spacecraft and venus states
% from spacecraft launch to encounter. It also compute the error between
% s/c and venus position and velocity, the date of arrival and the time of
% flight.
%                                                                          
% INPUT:
% -N_iter[1]: Maximum number of iteration for the while loop to find a
%             solution.
% - data[struct]: Struct containing adimentionalization constants and other
%                 spacecraft and environment parameters.
% -x0[6,1]: Spacecraft position and velocity at the time of launch.
% -varargin{1}: If present, it contains a non-random guess of the initial
%               costates of the spacecraft and the final time.
%
% OUTPUT:
% - sol[8,1]: Vector containing the computed initial costates and the final
%             time.
% - errors[struct]: Struct containing two scalar values representing the
%                   norm of the distance between venus and the s/c and the
%                   norm of the relative velocity at the final time.
% - T[n]: Vector containing the timestamps of the integration steps of the
%         s/c dynamics.
% -xx[n,14]: Matrix containing the 14 states and costates of the spacecraft
%            during the whole flight. The number of rows depends on 
%            integration options. 
% - venus[6,n]: Matrix containing venus position and velocity for each
%               column. The number of columns is not fix.
% -Timedata[struct]: struct containing a scalar value representing the time
%                    of flight in days and a char variable with the arrival
%                    date.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize exitflag and iteration counter
exitflag = 0;
i = 0;

% define initial guess and maximum thrust value in case of numerical
% continuation exploitation
if nargin > 3
    sol = varargin{1};
    data.Tmax = data.Tmax_vect;
end

% Iterate the initial guess untill a solution is found or the maximum
% number of iterations if reached
while exitflag ~=1 && i<N_iter

    % choose initial guess for both the standard case and the one which
    % exploits numerical continuation
    if nargin == 3
        guess = [40*rand(6,1) - 20;20*rand(1); 2*pi*rand(1)+data.t0];
    else
        guess = sol;
    end

    % solve zero-finding problem
    [sol,~,exitflag] = fsolve(@optimization,guess,[],x0,data);

    % increase iteration counter
    i = i+1;
end

% check if a proper solution was found
if i == N_iter && ~exitflag
    error('Algorithm did not converge')
end

% propagate s/c orbit with the correct initial costate and final time and
% compute venus states at the same time instants.
[T,xx] = ode113(@TPBVP,[data.t0,sol(8)],[x0;sol(1:7)],data.ode_opt,data);
venus = cspice_spkezr('Venus',T'*data.t,data.frame,'NONE',data.center);

% required time informations
tof = sol(8)-data.t0;
Timedata.tof = days(seconds(tof*data.t));
Timedata.tf = cspice_et2utc(sol(8)*data.t,'C',3);

% required errors
errors.pos = norm(venus(1:3,end)-xx(end,1:3)'*data.l);
errors.vel = norm(venus(4:6,end)-xx(end,4:6)'*data.l/data.t)*1e3;

end

function plotOrbits(venus,data,xx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the orbits of venus and the spacecraft.
%
% INPUT:
% - venus[6,n]: Matrix containing venus position and velocity for each
%               column. The number of columns is not fixed.
% - data[struct]: Struct containing adimentionalization constants and other
%                 spacecraft and environment parameters.
% -xx[n,14]: Matrix containing the 14 states and costates of the spacecraft
%            during the whole flight. The number of rows depends on 
%            integration options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % plot s/c trajectory
 figure()
 hold on
 grid on
 plot3(venus(1,end)/data.l,venus(2,end)/data.l,venus(3,end)/data.l,'diamond','MarkerFaceColor',...
     "#D95319",'MarkerSize',15,'DisplayName','Venus at encounter')
 plot3(venus(1,1)/data.l,venus(2,1)/data.l,venus(3,1)/data.l,'square','MarkerFaceColor',...
     "#D95319",'MarkerSize',15,'DisplayName','Venus at s/c departure')
 plot3(venus(1,:)/data.l,venus(2,:)/data.l,venus(3,:)/data.l,'Color',...
     "#D95319",'DisplayName','Venus orbit')
 plot3(xx(:,1),xx(:,2),xx(:,3),'Color', "#0072BD",'DisplayName','S/C trajectory')
 plot3(xx(end,1),xx(end,2),xx(end,3),'o','MarkerFaceColor',...
     "#0072BD",'MarkerSize',10,'DisplayName','S/C at encounter')
 plot3(0,0,0,'o','MarkerFaceColor',...
     "#EDB120",'MarkerEdgeColor','none','MarkerSize',15,'DisplayName','Sun')

 legend('Location','best','Orientation','vertical');
 xlabel('x [AU]')
 ylabel('y [AU]')
 zlabel('z [AU]')
view(30,20)
 axis equal
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
