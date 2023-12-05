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

%% Set initial conditions

t_sep = cspice_str2et('2010-08-12T05:27:39.114');
r1_0  = [4622.232026629; 5399.3369588058; -0.0212138165769957];
v1_0  = [0.812221125483763; -0.721512914578826; 7.42665302729053]; 
r2_0  = [4621.69343340281; 5399.26386352847; -3.09039248714313];
v2_0  = [0.813960847513811; -0.719449862738607; 7.42706066911294];
P0 = [ +5.6e-7, +3.5e-7, -7.1e-8,    0,        0,        0;
       +3.5e-7, +9.7e-7, +7.6e-8,    0,        0,        0;
       -7.1e-8, +7.6e-8, +8.1e-8,    0,        0,        0;
           0,       0,       0,  +2.8e-11,     0,        0;
           0,       0,       0,      0,    +2.7e-11,     0;
           0,       0,       0,      0,        0,    +9.6e-12 ];
phi0 = reshape(eye(6),36,1);

%% Point 1

% compute orbital period of satellite 1
settings.mu = cspice_bodvrd('Earth','GM',1);
a = -0.5 * settings.mu / ( norm(v1_0)^2/2 - settings.mu/norm(r1_0));
settings.T1 = 2*pi*sqrt(a^3/settings.mu);
N = 10;
n = 6;

settings.ut.alpha = 0.1;
settings.ut.beta = 2;

Mango_lin.states = [ [r1_0;v1_0;phi0] , zeros(42,N)];
Mango_lin.P = [zeros(n,n,N+1)];
Mango_lin.P(:,:,1) = P0;
Tango_lin.states = [ [r2_0;v2_0;phi0] , zeros(42,N)];
Tango_lin.P = Mango_lin.P;
settings.ode_opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

Mango_ut.states = [ [r1_0;v1_0;] , zeros(6,N)];
Mango_ut.P = Mango_lin.P;
Mango_ut.sigma_points = false;
Tango_ut.states = [ [r2_0;v2_0;] , zeros(6,N)];
Tango_ut.P = Tango_lin.P;
Tango_ut.sigma_points = false;

for i = 2:N+1
    % LinCov
    [Mango_lin.states(:,i),Mango_lin.P(:,:,i)] = ...
        LinCov([Mango_lin.states(1:6,i-1);phi0], Mango_lin.P(:,:,i-1), settings);
    [Tango_lin.states(:,i),Tango_lin.P(:,:,i)] = ...
        LinCov([Tango_lin.states(1:6,i-1);phi0], Tango_lin.P(:,:,i-1), settings);

    % UT
    [Mango_ut.states(:,i),Mango_ut.P(:,:,i),Mango_ut.sigma_points] = ...
        UT(Mango_ut.states(:,i-1), Mango_ut.P(:,:,i-1),Mango_ut.sigma_points, settings);
    [Tango_ut.states(:,i),Tango_ut.P(:,:,i),Tango_ut.sigma_points] = ...
        UT(Tango_ut.states(:,i-1), Tango_ut.P(:,:,i-1),Mango_ut.sigma_points, settings);
end

%% Ex 2
delta_r_lin = zeros(N+1,1);
P_sum_lin = zeros(6,6,N+1);
delta_r_lin_lim =  zeros(N+1,1);

delta_r_ut= zeros(N+1,1);
P_sum_ut = zeros(6,6,N+1);
delta_r_ut_lim =  zeros(N+1,1);

flag.found_delta_r_lin = false;
flag.found_delta_r_ut = false;

for i = 1:N+1
    delta_r_lin(i) = norm( Mango_lin.states(1:3,i) - Tango_lin.states(1:3,i) );
    delta_r_ut(i) = norm( Mango_ut.states(1:3,i) - Tango_ut.states(1:3,i) );
    
    P_sum_lin(:,:,i) = Tango_lin.P(:,:,i) + Mango_lin.P(:,:,i);
    P_sum_ut(:,:,i) = Tango_ut.P(:,:,i) + Mango_ut.P(:,:,i);
    
    delta_r_lin_lim(i) = 3 * sqrt( max( eig( P_sum_lin(:,:,i) ) ) );
    delta_r_ut_lim(i) = 3 * sqrt( max( eig( P_sum_ut(:,:,i) ) ) );

    if delta_r_lin(i) < delta_r_lin_lim(i) && ~flag.found_delta_r_lin
        flag.found_delta_r_lin = true;
        N_lin_lim = i-1;
    end

    if delta_r_ut(i) < delta_r_ut_lim(i) && ~flag.found_delta_r_ut
        flag.found_delta_r_ut = true;
        N_lin_ut = i-1;
    end
end


%% functions
function [dxx] = TBP(~,xx,mu)

x = xx(1);
y = xx(2);
z = xx(3);
r_norm = norm([x;y;z]);

if length(xx) == 42
    % Put the STM PHI in matrix form
    phi = reshape(xx(7:end),6,6);

    % Jacobian of f
    A = zeros(6);
    A(1:3,4:6) = eye(3);
    A(4,1) = -(mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
    A(4,2) = (3*mu*x*y)/(r_norm)^5;
    A(4,3) = (3*mu*x*z)/(r_norm)^5;
    A(5,1) = (3*mu*x*y)/(r_norm)^5;
    A(5,2) = -(mu*(x^2 - 2*y^2 + z^2))/(r_norm)^5;
    A(5,3) = (3*mu*y*z)/r_norm^5;
    A(6,1) = (3*mu*x*z)/r_norm^5;
    A(6,2) = (3*mu*y*z)/r_norm^5;
    A(6,3) = -(mu*(x^2 + y^2 - 2*z^2))/r_norm^5;

    % Compute the derivative of the STM
    dphi = A*phi;
    dxx(7:42) = reshape(dphi,36,1);
end

% RHS
dxx(1:3) = xx(4:6);
dxx(4:6) = -mu/r_norm^3 * [x;y;z];

dxx = dxx';

end

function [x,P] = LinCov(x0,P0,settings)

mu = settings.mu;
ode_opt = settings.ode_opt;
tf = settings.T1;

[~,x] = ode78(@TBP,[0,tf],x0,ode_opt,mu);

x = x(end,:);
phi = reshape(x(7:end),6,6);
P = phi*P0*phi';
end

function [y_mean,P,sigma_points] = UT(x0,P0,sigma_points,settings)

ode_opt = settings.ode_opt;
tf = settings.T1;
mu = settings.mu;
alpha = settings.ut.alpha;
beta = settings.ut.beta;

n = length(x0);
lambda = alpha^2*n-n;
mat = sqrtm((n+lambda)*P0);

gamma = zeros(n,2*n+1);
weight_mean = zeros(2*n+1,1);
weight_cov = zeros(2*n+1,1);
y_mean = zeros(n,1);
P = zeros(6);

if ~sigma_points
chi = zeros(n,2*n+1);
chi(:,1) = x0;
for i = 1:n
    chi(:,i+1) = x0 + mat(:,i);
    chi(:,i+1+n) = x0 - mat(:,i);
end
else
    chi=sigma_points;
end

weight_mean(1) = lambda/(n+lambda);
weight_cov(1) = weight_mean(1) + (1-alpha^2+beta);

weight_mean(2:end) = 1/(2*(n+lambda)) * ones(2*n,1);
weight_cov(2:end) = weight_mean(2:end);

for i = 1:2*n+1
    [~,x] = ode78(@TBP,[0,tf],chi(:,i),ode_opt,mu);
    gamma(:,i) = x(end,:)';
    y_mean = y_mean + weight_mean(i)*gamma(:,i);
end

for i = 1:2*n+1
    P = P + weight_cov(i)*(gamma(:,i)-y_mean)*(gamma(:,i)-y_mean)';
end

sigma_points = gamma;

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