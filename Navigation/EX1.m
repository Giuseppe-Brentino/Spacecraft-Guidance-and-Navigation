% Spacecraft Guidance and Navigation (2023/2024)
% Assignment #2: Exercise 1
% Author: Giuseppe Brentino

clearvars; close all; clc;
addpath('.\kernels\')
cspice_furnsh('assignment02.tm')
plotStyle;

rng default
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

% Set up initial parameters for Earth
settings.mu = cspice_bodvrd('Earth', 'GM', 1);

% Compute the semi-major axis
a = -0.5 * settings.mu / (norm(v1_0)^2 / 2 - settings.mu / norm(r1_0));

% Compute the orbital period 'T1' of satellite 1
settings.T1 = 2 * pi * sqrt(a^3 / settings.mu);

% Number of iterations and dimensions
N = 10;
n = 6;

% Unscented Transform parameters
settings.ut.alpha = 0.1;
settings.ut.beta = 2;

% Initial state vectors and covariance matrices for Mango and Tango satellites
Mango_lin.states = [[r1_0; v1_0; phi0], zeros(42, N)];
Mango_lin.P = zeros(n, n, N + 1);
Mango_lin.P(:, :, 1) = P0;
Tango_lin.states = [[r2_0; v2_0; phi0], zeros(42, N)];
Tango_lin.P = Mango_lin.P;

% Options for ODE solver
settings.ode_opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-13);

% Unscented Transform states and covariance for Mango and Tango
Mango_ut.states = [[r1_0; v1_0;], zeros(6, N)];
Mango_ut.P = Mango_lin.P;
Mango_ut.sigma_points = [];
Tango_ut.states = [[r2_0; v2_0;], zeros(6, N)];
Tango_ut.P = Tango_lin.P;
Tango_ut.sigma_points = [];

% Iteratively update states and covariances using LinCov and Unscented Transform
for i = 2:N + 1
    % LinCov
    [Mango_lin.states(:, i), Mango_lin.P(:, :, i)] = ...
        LinCov([Mango_lin.states(1:6, i - 1); phi0], Mango_lin.P(:, :, i - 1), settings);
    [Tango_lin.states(:, i), Tango_lin.P(:, :, i)] = ...
        LinCov([Tango_lin.states(1:6, i - 1); phi0], Tango_lin.P(:, :, i - 1), settings);

    % Unscented Transform
    [Mango_ut.states(:, i), Mango_ut.P(:, :, i), Mango_ut.sigma_points] = ...
        UT(Mango_ut.states(:, i - 1), Mango_ut.P(:, :, i - 1), Mango_ut.sigma_points, settings);
    [Tango_ut.states(:, i), Tango_ut.P(:, :, i), Tango_ut.sigma_points] = ...
        UT(Tango_ut.states(:, i - 1), Tango_ut.P(:, :, i - 1), Tango_ut.sigma_points, settings);
end

%% Ex 2

% Initialize arrays and flags
delta_r_lin = zeros(N + 1, 1);
P_sum_lin = zeros(3, 3, N + 1);
delta_r_lin_lim = zeros(N + 1, 1);

delta_r_ut = zeros(N + 1, 1);
P_sum_ut = zeros(3, 3, N + 1);
delta_r_ut_lim = zeros(N + 1, 1);

flag.found_delta_r_lin = false;
flag.found_delta_r_ut = false;

% Loop over each orbit
for i = 1:N + 1
    % Calculate differences in position vectors
    delta_r_lin(i) = norm(Mango_lin.states(1:3, i) - Tango_lin.states(1:3, i));
    delta_r_ut(i) = norm(Mango_ut.states(1:3, i) - Tango_ut.states(1:3, i));

    % Sum of covariance matrices for position vectors
    P_sum_lin(:, :, i) = Tango_lin.P(1:3, 1:3, i) + Mango_lin.P(1:3, 1:3, i);
    P_sum_ut(:, :, i) = Tango_ut.P(1:3, 1:3, i) + Mango_ut.P(1:3, 1:3, i);

    % Limit for linearized error
    delta_r_lin_lim(i) = 3 * sqrt(eigs(P_sum_lin(:, :, i), 1));

    % Limit for Unscented Transform error
    delta_r_ut_lim(i) = 3 * sqrt(eigs(P_sum_ut(:, :, i), 1));

    % Check if linearized error is within limit
    if delta_r_lin(i) < delta_r_lin_lim(i) && ~flag.found_delta_r_lin
        flag.found_delta_r_lin = true;
        N_lin_lim = i - 1;
    end

    % Check if Unscented Transform error is within limit
    if delta_r_ut(i) < delta_r_ut_lim(i) && ~flag.found_delta_r_ut
        flag.found_delta_r_ut = true;
        N_ut_lim = i - 1;
    end
end

figure()
hold on
grid on
plot(0:10,delta_r_lin)
plot(0:10,delta_r_lin_lim)
xline(N_lin_lim)
legend('$\Delta$r','$\Delta$r lim', 'First unsafe orbit')
xlabel('Number of orbits [-]')
ylabel('$\Delta$r [km]')


figure()
hold on
grid on
plot(0:10,delta_r_ut)
plot(0:10,delta_r_ut_lim)
xline(N_ut_lim)
legend('$\Delta$r','$\Delta$r lim', 'First unsafe orbit')
xlabel('Number of orbits [-]')
ylabel('$\Delta$r [km]')
%% Ex 3

% Set the number of simulations
settings.n_sim = 500;

% Generate random initial conditions for Mango
Mango_mc.x_0 = mvnrnd([r1_0; v1_0], P0, settings.n_sim)';

% Generate random initial conditions for Tango
Tango_mc.x_0 = mvnrnd([r2_0; v2_0], P0, settings.n_sim)';

% Initialize mean and covariance matrices for Mango
Mango_mc.x_mean = [[r1_0; v1_0], zeros(6, N)];
Mango_mc.P = zeros(6, 6, N + 1);
Mango_mc.P(:, :, 1) = P0;

% Initialize mean and covariance matrices for Tango
Tango_mc.x_mean = [[r2_0; v2_0], zeros(6, N)];
Tango_mc.P = Mango_mc.P;

% Monte Carlo simulations loop
for i = 2:N + 1
    % Perform Monte Carlo simulation for Mango
    [Mango_mc.x_mean(:, i), Mango_mc.P(:, :, i), states] = MonteCarlo(Mango_mc.x_0, settings);
    Mango_mc.x_0 = states;

    % Perform Monte Carlo simulation for Tango
    [Tango_mc.x_mean(:, i), Tango_mc.P(:, :, i), states] = MonteCarlo(Tango_mc.x_0, settings);
    Tango_mc.x_0 = states;
end

%%% Plots

% Covariance time evolution

Mango_lin.Pr_eig = zeros(N+1,1);
Mango_lin.Pv_eig = zeros(N+1,1);
Mango_ut.Pr_eig = zeros(N+1,1);
Mango_ut.Pv_eig = zeros(N+1,1);
Mango_mc.Pr_eig = zeros(N+1,1);
Mango_mc.Pv_eig = zeros(N+1,1);

Tango_lin.Pr_eig = zeros(N+1,1);
Tango_lin.Pv_eig = zeros(N+1,1);
Tango_ut.Pr_eig = zeros(N+1,1);
Tango_ut.Pv_eig = zeros(N+1,1);
Tango_mc.Pr_eig = zeros(N+1,1);
Tango_mc.Pv_eig = zeros(N+1,1);

for i = 1:N+1    
    Mango_mc.Pr_eig(i) = 3*sqrt( eigs(Mango_mc.P(1:3,1:3,i),1) );
    Mango_mc.Pv_eig(i) = 3*sqrt( eigs(Mango_mc.P(4:6,4:6,i),1) );
    Mango_ut.Pr_eig(i) = 3*sqrt( eigs(Mango_ut.P(1:3,1:3,i),1) );
    Mango_ut.Pv_eig(i) = 3*sqrt( eigs(Mango_ut.P(4:6,4:6,i),1) );
    Mango_lin.Pr_eig(i) = 3*sqrt( eigs(Mango_lin.P(1:3,1:3,i),1) );
    Mango_lin.Pv_eig(i) = 3*sqrt( eigs(Mango_lin.P(4:6,4:6,i),1) );

    Tango_mc.Pr_eig(i) = 3*sqrt( eigs(Tango_mc.P(1:3,1:3,i),1) );
    Tango_mc.Pv_eig(i) = 3*sqrt( eigs(Tango_mc.P(4:6,4:6,i),1) );
    Tango_ut.Pr_eig(i) = 3*sqrt( eigs(Tango_ut.P(1:3,1:3,i),1) );
    Tango_ut.Pv_eig(i) = 3*sqrt( eigs(Tango_ut.P(4:6,4:6,i),1) );
    Tango_lin.Pr_eig(i) = 3*sqrt( eigs(Tango_lin.P(1:3,1:3,i),1) );
    Tango_lin.Pv_eig(i) = 3*sqrt( eigs(Tango_lin.P(4:6,4:6,i),1) );
end

orbit_index = 0:N;
figure()
hold on
grid on
subplot(2,2,1)
hold on
plot(orbit_index, Mango_lin.Pr_eig)
plot(orbit_index, Mango_ut.Pr_eig,'--')
plot(orbit_index, Mango_mc.Pr_eig,'-.')
title('Mango Pr')
legend('LinCov','UT','MC')
xlabel('Number of orbits [-]')
ylabel('$\sigma_r$')
subplot(2,2,2)
hold on
plot(orbit_index, Mango_lin.Pv_eig)
plot(orbit_index, Mango_ut.Pv_eig,'--')
plot(orbit_index, Mango_mc.Pv_eig,'-.')
title('Mango Pv')
legend('LinCov','UT','MC')
xlabel('Number of orbits [-]')
ylabel('$\sigma_v$')
subplot(2,2,3)
hold on
plot(orbit_index, Tango_lin.Pr_eig)
plot(orbit_index, Tango_ut.Pr_eig,'--')
plot(orbit_index, Tango_mc.Pr_eig,'-.')
title('Tango Pr')
legend('LinCov','UT','MC')
xlabel('Number of orbits [-]')
ylabel('$\sigma_r$')
subplot(2,2,4)
hold on
plot(orbit_index, Tango_lin.Pv_eig)
plot(orbit_index, Tango_ut.Pv_eig,'--')
plot(orbit_index, Tango_mc.Pv_eig,'-.')
xlabel('Number of orbits [-]')
ylabel('$\sigma_v$')
title('Tango Pv')
legend('LinCov','UT','MC')

% Covariance ellipses Mango
i_v = r1_0/norm(r1_0);
k_v = cross(r1_0,v1_0)/norm(cross(r1_0,v1_0));
j_v = cross(k_v,i_v);
R = [i_v';j_v';k_v'];
Mango_mc.rotated_pos = zeros(3,settings.n_sim);

for i = 1:settings.n_sim
    Mango_mc.rotated_pos(:,i) = R * Mango_mc.x_0(1:3,i);
end

Mango_mc.rotated_mean = R * Mango_mc.x_mean(1:3,end);
Mango_lin.rotated_mean = R * Mango_lin.states(1:3,end);
Mango_ut.rotated_mean = R * Mango_ut.states(1:3,end);

Mango_mc.rotated_cov = R*Mango_mc.P(1:3,1:3,end)*R';
Mango_lin.rotated_cov = R*Mango_lin.P(1:3,1:3,end)*R';
Mango_ut.rotated_cov = R*Mango_ut.P(1:3,1:3,end)*R';

figure()
hold on
grid on
plot(Mango_mc.rotated_pos(1,:),Mango_mc.rotated_pos(2,:),'k.');
drawEllipse(Mango_lin.rotated_mean(1:2),Mango_lin.rotated_cov(1:2,1:2),3);
drawEllipse(Mango_ut.rotated_mean(1:2),Mango_ut.rotated_cov(1:2,1:2),3);
drawEllipse(Mango_mc.rotated_mean(1:2),Mango_mc.rotated_cov(1:2,1:2),3);
plot(Mango_lin.rotated_mean(1),Mango_lin.rotated_mean(2),'o','MarkerSize',10,'MarkerFaceColor','none')
plot(Mango_ut.rotated_mean(1),Mango_ut.rotated_mean(2),'*','MarkerSize',10,'MarkerFaceColor','none')
plot(Mango_mc.rotated_mean(1),Mango_mc.rotated_mean(2),'s','MarkerSize',10,'MarkerFaceColor','none')
xlabel('x [km]')
ylabel('y [km]')
legend('Propagated points','LinCov ellipse','UT ellipse','MC ellipse',...
    'LinCov mean','UT mean','MC mean')

% Covariance ellipses Tango
i_v = r2_0/norm(r2_0);
k_v = cross(r2_0,v2_0)/norm(cross(r2_0,v2_0));
j_v = cross(k_v,i_v);
R = [i_v';j_v';k_v'];
Tango_mc.rotated_pos = zeros(3,settings.n_sim);

for i = 1:settings.n_sim
    Tango_mc.rotated_pos(:,i) = R * Tango_mc.x_0(1:3,i);
end

Tango_mc.rotated_mean = R * Tango_mc.x_mean(1:3,end);
Tango_lin.rotated_mean = R * Tango_lin.states(1:3,end);
Tango_ut.rotated_mean = R * Tango_ut.states(1:3,end);

Tango_mc.rotated_cov = R*Tango_mc.P(1:3,1:3,end)*R';
Tango_lin.rotated_cov = R*Tango_lin.P(1:3,1:3,end)*R';
Tango_ut.rotated_cov = R*Tango_ut.P(1:3,1:3,end)*R';

figure()
hold on
grid on
plot(Tango_mc.rotated_pos(1,:),Tango_mc.rotated_pos(2,:),'k.');
drawEllipse(Tango_lin.rotated_mean(1:2,end),Tango_lin.rotated_cov(1:2,1:2),3);
drawEllipse(Tango_ut.rotated_mean(1:2,end),Tango_ut.rotated_cov(1:2,1:2),3);
drawEllipse(Tango_mc.rotated_mean(1:2,end),Tango_mc.rotated_cov(1:2,1:2),3);
plot(Tango_lin.rotated_mean(1),Tango_lin.rotated_mean(2),'o','MarkerSize',10,'MarkerFaceColor','none')
plot(Tango_ut.rotated_mean(1),Tango_ut.rotated_mean(2),'*','MarkerSize',10,'MarkerFaceColor','none')
plot(Tango_mc.rotated_mean(1),Tango_mc.rotated_mean(2),'s','MarkerSize',10,'MarkerFaceColor','none')
xlabel('x [km]')
ylabel('y [km]')
legend('Propagated points','LinCov ellipse','UT ellipse','MC ellipse',...
    'LinCov mean','UT mean','MC mean')

%% functions

function [dxx] = TBP(~, xx, mu)
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
r_norm = norm([x; y; z]);

% Check if the state vector includes the STM
if length(xx) == 42
    % Extract the STM PHI from the state vector
    phi = reshape(xx(7:end), 6, 6);

    % Jacobian matrix of the TBP equations
    A = zeros(6);
    A(1:3, 4:6) = eye(3);
    A(4, 1) = -(mu * (-2 * x^2 + y^2 + z^2)) / (r_norm^5);
    A(4, 2) = (3 * mu * x * y) / r_norm^5;
    A(4, 3) = (3 * mu * x * z) / r_norm^5;
    A(5, 1) = (3 * mu * x * y) / r_norm^5;
    A(5, 2) = -(mu * (x^2 - 2 * y^2 + z^2)) / r_norm^5;
    A(5, 3) = (3 * mu * y * z) / r_norm^5;
    A(6, 1) = (3 * mu * x * z) / r_norm^5;
    A(6, 2) = (3 * mu * y * z) / r_norm^5;
    A(6, 3) = -(mu * (x^2 + y^2 - 2 * z^2)) / r_norm^5;

    % Compute the derivative of the STM
    dphi = A * phi;
    dxx(7:42) = reshape(dphi, 36, 1);
end

% Right-hand side (RHS) of the TBP equations
dxx(1:3) = xx(4:6);
dxx(4:6) = -mu / r_norm^3 * [x; y; z];

% Transpose to ensure a column vector
dxx = dxx';

end

function [x,P] = LinCov(x0,P0,settings)
% LinCov - Propagate mean state and covariance matrix through a linear 
%          transformation.
%
%   Inputs:
%       - x0: Initial state vector.
%             Column vector.
%       - P0: Initial covariance matrix.
%             Square matrix.
%       - settings: Struct containing the following parameters:
%           - settings.mu: Gravitational parameter.
%           - settings.ode_opt: Options for the ODE solver.
%           - settings.T1: Final time for integration.
%
%   Outputs:
%       - x: Final state vector after integration.
%            Column vector.
%       - P: Final covariance matrix.
%            Square matrix.
%

mu = settings.mu;
ode_opt = settings.ode_opt;
tf = settings.T1;

% Solve the Two-Body Problem
[~, x] = ode78(@TBP, [0, tf], x0, ode_opt, mu);

% Extract the final state vector
x = x(end, :);

% Extract the STM from the state vector
phi = reshape(x(7:end), 6, 6);

% Compute the final covariance matrix using the linearized STM
P = phi * P0 * phi';

end

function [y_mean,P,sigma_points] = UT(x0,P0,sigma_points,settings)
% UT - Unscented Transform for estimating the mean and covariance of the
%      state vector using the Two-Body Problem (TBP) model.
%
%   Inputs:
%       - x0: Initial state vector.
%             Column vector.
%       - P0: Initial covariance matrix.
%             Square matrix.
%       - sigma_points: Sigma points for the Unscented Transform (optional,
%                       it can also be an empty matrix).
%                        Matrix with n rows and 2n+1 columns, where n is the
%                        dimension of the state vector.
%       - settings: Struct containing the following parameters:
%           - settings.ode_opt: Options for the ODE solver.
%           - settings.T1: Final time for integration.
%           - settings.mu: Gravitational parameter.
%           - settings.ut.alpha: Scaling parameter for sigma points.
%           - settings.ut.beta: Incorporates prior knowledge about the
%                                distribution of the state.
%
%   Outputs:
%       - y_mean: Estimated mean of the state vector.
%                 Column vector.
%       - P: Estimated covariance matrix of the state vector.
%            Square matrix.
%       - sigma_points: Updated or computed sigma points.
%                       Matrix with n rows and 2n+1 columns.
%

ode_opt = settings.ode_opt;
tf = settings.T1;
mu = settings.mu;
alpha = settings.ut.alpha;
beta = settings.ut.beta;

n = length(x0);
lambda = alpha^2 * n - n;
mat = sqrtm((n + lambda) * P0);

weight_mean = zeros(2 * n + 1, 1);
weight_cov = zeros(2 * n + 1, 1);
y_mean = zeros(n, 1);
P = zeros(n);

% Check if sigma_points are provided; if not, generate them
if isempty(sigma_points)
    chi = zeros(n, 2 * n + 1);
    chi(:, 1) = x0;
    for i = 1:n
        chi(:, i + 1) = x0 + mat(:, i);
        chi(:, i + 1 + n) = x0 - mat(:, i);
    end
else
    chi = sigma_points;
end

% Compute weights for mean and covariance
weight_mean(1) = lambda / (n + lambda);
weight_cov(1) = weight_mean(1) + (1 - alpha^2 + beta);

weight_mean(2:end) = 1 / (2 * (n + lambda)) * ones(2 * n, 1);
weight_cov(2:end) = weight_mean(2:end);

% Propagate sigma points through the dynamics model
for i = 1:2 * n + 1
    [~, x] = ode78(@TBP, [0, tf], chi(:, i), ode_opt, mu);
    sigma_points(:, i) = x(end, :)';
    y_mean = y_mean + weight_mean(i) * sigma_points(:, i);
end

% Compute the covariance matrix P
for i = 1:2 * n + 1
    P = P + weight_cov(i) * (sigma_points(:, i) - y_mean) * (sigma_points(:, i) - y_mean)';
end


end

function [x_mean, covariance, states] = MonteCarlo(x0,settings)
% MonteCarlo - Performs Monte Carlo simulation for the Two-Body Problem (TBP).
%
%   Inputs:
%       - x0: Initial state vectors for each simulation.
%             Matrix with 6 rows and n_sim columns.
%       - settings: Struct containing the following parameters:
%           - settings.ode_opt: Options for the ODE solver.
%           - settings.T1: Final time for integration.
%           - settings.mu: Gravitational parameter.
%           - settings.n_sim: Number of simulations to run.
%
%   Outputs:
%       - x_mean: Mean state vector over all simulations.
%                 Column vector.
%       - covariance: Covariance matrix of the state vectors.
%                     Square matrix.
%       - states: Matrix containing the final state vectors for each simulation.
%                 Each column represents a simulation, and each row represents
%                 a component of the state vector.
%

ode_opt = settings.ode_opt;
tf = settings.T1;
mu = settings.mu;
n_sim = settings.n_sim;

% Preallocate matrix to store final state vectors for each simulation
states = zeros(6, n_sim);

% Parallel loop for running simulations
parfor i = 1:n_sim
    [~, x] = ode78(@TBP, [0, tf], x0(:, i), ode_opt, mu);
    states(:, i) = x(end, :)';
end

% Calculate mean and covariance over all simulations
x_mean = mean(states, 2);
covariance = cov(states');

end

function [varargout] = drawEllipse(mean,P,n_sigma)
% drawEllipse - Draws an ellipse based on mean and covariance matrix.
%
%   Inputs:
%       - mean: Mean vector of the distribution.
%               Column vector.
%       - P: Covariance matrix.
%            Square matrix.
%       - n_sigma: Number of standard deviations for ellipse size.
%                  Positive scalar.
%
%   Outputs:
%       - If called with no output arguments, the function plots the ellipse.
%       - If called with two output arguments, the function returns the x and y
%         coordinates of the ellipse without plotting.
%

% Build circle
n_points = 200;
alpha = 2 * pi / n_points * (0:n_points);
circle = [cos(alpha); sin(alpha)];

% Singular Value Decomposition (SVD) to find ellipse parameters
[R, D] = svd(P);
d = sqrt(D);

% Transform the circle to an aligned ellipse, then rotate it, and finally
% scale it based on the specified number of standard deviations
ellipse = n_sigma * R * d * circle;

% Shift the ellipse to the mean position
x = mean(1) + ellipse(1, :);
y = mean(2) + ellipse(2, :);

if nargout == 2
    varargout{1} = x;
    varargout{2} = y;
else
    plot(x, y);
end

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
set(0,'defaultLineMarkerFaceColor','auto')
set(0,'defaultLineMarkerEdgeColor','auto')
% legend:
set(0, 'defaultLegendLocation','southoutside');
set(0, 'defaultLegendOrientation','horizontal');
set(0, 'defaultLegendFontSize',12);
% axes:
set(0,'defaultAxesFontSize',16);
end