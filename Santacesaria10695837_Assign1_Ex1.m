% Spacecraft Guidance and Navigations (2023/2024)
% Assignment 1, Excercise 1
% Author: Matteo Santacesaria

% add path in order for cspice to work properly
run micepath.m;

%% Ex 2
clear; close all; clc;

% Set mu for the 3D Earth-Moon CR3BP 
mu = 0.012150;

% Initial conditions
x0  = 0.829380407710981;
y0  = 0;
z0  = 0.0901786827424052;
vx0 = 0;
vy0 = 0.181868545812666;
vz0 = 0;

% Vector of initial conditions
xx0 = [x0;y0;z0;vx0;vy0;vz0];

% Propagation time
tf = 100.0;

% flag set to false to continue the propagation up to final time
[~,~,xx_F]  = propagate(0,xx0,tf,mu,false); 
[~,~,xx_B]  = propagate(0,xx0,-tf,mu,false);

% Initialization
err_vxf = 1;    % Set this to a high value not to stop at first iteration
err_vzf = 1;
Nmax    = 150;  % Set a maximum number of iterations to avoid it to get stuck
iter    = 0;    % Initialize the iteration counter
tol     = 1e-8; % Set the desired tolerance

e = zeros(length(xx0));
PHI_N = zeros(length(xx0));

% Set the update to the first guess of the velocity
vy0_new = vy0;
z0_new  = z0;
vxf_ref = 0;
vzf_ref = 0;

while (abs(err_vxf)>tol || abs(err_vzf)>tol) && iter < Nmax
    xx0_new = [x0;y0;z0_new;vx0;vy0_new;vz0];
    [te,xf_nom]   = propagate(0,[x0;y0;z0_new;vx0;vy0_new;vz0],tf,mu);
    
    % Computation of the STM by columnbs
    for i = 1 : length(xx0_new)
        e(i,i) = max(sqrt(eps), sqrt(eps)*abs(xx0_new(i)));     % Perturbation
        [te,xf,~]  = propagate(0,xx0_new+e(:,i),tf,mu);         % Perturbed solution
        PHI_N(:,i) = (xf - xf_nom) / e(i,i);                    % Forward differences
    end

    % Compute the deviation in the final state (x velocity)
    err_vxf = abs(xf(4) - vxf_ref);
    err_vzf = abs(xf(6) - vzf_ref);
    PHI_vv = [PHI_N(4, 3), PHI_N(4, 5); PHI_N(6, 3), PHI_N(6, 5)];
    
    dv = PHI_vv\[-xf(4); -xf(6)];

    % Compute the correction
    vy0_new = dv(2) + vy0_new;
    z0_new = dv(1) + z0_new;
%     t_out=te;
    % Update iteration counter
    iter = iter+1;
end

xx0_new1 = [x0;y0;z0_new;vx0;vy0_new;vz0];
[te,xf] = propagate(0,xx0_new1,2*te,mu,false);
err_x0 = (xf-xx0_new);
fprintf('\nvx(t_e) = %.3e\nNumber of iterations: %d',xf(4),iter);
fprintf('\nvz(t_e) = %.3e\nNumber of iterations: %d',xf(6),iter);
[te, ~, xx_F] = propagate(0,xx0_new, te,mu,false);
[te, ~, xx_B] = propagate(0,xx0_new, -te,mu,false);

%% plots

figure(1)
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
hold on
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',2)
plot(1-mu,0,'o','Marker', 'o', 'MarkerFaceColor','y','MarkerSize',8)
plot(-mu,0,'o','Marker', 'o', 'MarkerFaceColor','g','MarkerSize',8)
xlabel('x [-]','FontSize',14)
ylabel('y [-]','FontSize',14)
zlabel('z[-]','FontSize',14)
legend(["Forward","Backward"],'Location','best','FontSize',14)

figure(2)
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),LineWidth=2)
hold on
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),LineWidth=2)
plot(1-mu,0,'o','Marker', 'o', 'MarkerFaceColor','y','MarkerSize',8)
plot(-mu,0,'o','Marker', 'o', 'MarkerFaceColor','g','MarkerSize',8)
xlabel('x [-]','FontSize',14)
ylabel('y [-]','FontSize',14)
zlabel('z[-]','FontSize',14)
legend(["Forward","Backward"],'Location','best','FontSize',14)

%% Ex 3
clear; close all; clc;
clc, clearvars, close all
n=100;
x0=linspace(0.829380407710981,0.873,n);
y0  = 0;
z0  = 0.0901786827424052;
vx0 = 0;
vy0 = 0.181868545812666;
vz0 = 0;
mu = 0.012150; %
tf = 2.0;   % Propagation time
for j=1:length(x0)
    xx0 = [x0(j);y0;z0;vx0;vy0;vz0];

%     figure(1)
    [~,~,xx_F]  = propagate(0,xx0,tf,mu,false);  % flag set to false to continue the propagation up to final time
    [~,~,xx_B]  = propagate(0,xx0,-tf,mu,false); % flag set to false to continue the propagation up to final time
    figure(3)
    plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
    hold on
    plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',2)
    plot(1-mu,0,'o','Marker', 'o', 'MarkerFaceColor','y','MarkerSize',8)
    plot(-mu,0,'o','Marker', 'o', 'MarkerFaceColor','g','MarkerSize',8)
    xlabel('x [-]','FontSize',14)
    ylabel('y [-]','FontSize',14)
    zlabel('z[-]','FontSize',14)
%     legend(["Forward","Backward"],'Location','best',FontSize=14)
    % Initialization
    err_vxf = 1; % Set this to a high value not to stop at first iteration
    err_vzf = 1;
    Nmax    = 150;   % Set a maximum number of iterations to avoid it to get stuck
    iter    = 0;    % Initialize the iteration counter
    tol     = 1e-9; % Set the desired tolerance
    e = zeros(length(xx0));
    PHI_N = zeros(length(xx0));

    % Set the update to the first guess of the velocity
    vy0_new = vy0;
    z0_new  = z0;
    vxf_ref = 0;
    vzf_ref = 0;
    while (abs(err_vxf)>tol || abs(err_vzf)>tol) && iter < Nmax
        xx0_new = [x0(j);y0;z0_new;vx0;vy0_new;vz0];
        [te,xf_nom]   = propagate(0,[x0(j);y0;z0_new;vx0;vy0_new;vz0],tf,mu);

        %Computation of the STM by columnbs
        for i = 1 : length(xx0_new)
            e(i,i) = max(sqrt(eps), sqrt(eps)*abs(xx0_new(i)));       % Perturbation
            [te,xf,~]  = propagate(0,xx0_new+e(:,i),tf,mu);          % Perturbed solution
            PHI_N(:,i) = (xf - xf_nom) / e(i,i);               % Forward differences
        end
        % Compute the deviation in the final state (x velocity)
        err_vxf = abs(xf(4) - vxf_ref);
        err_vzf = abs(xf(6) - vzf_ref);
        PHI_vv = [PHI_N(4, 3), PHI_N(4, 5); PHI_N(6, 3), PHI_N(6, 5)];

        dv = PHI_vv\[-xf(4); -xf(6)];

        % Compute the correction
        vy0_new = dv(2) + vy0_new;
        z0_new = dv(1) + z0_new;
        % Update iteration counter
        iter = iter+1;
    end

    xx0_new1 = [x0(j);y0;z0_new;vx0;vy0_new;vz0];
    [te,xf] = propagate(0,xx0_new1,2*te,mu,false);
    err_x0 = (xf-xx0_new);
    fprintf('\nvx(t_e) = %.3e\nNumber of iterations: %d',xf(4),iter);
    fprintf('\nvz(t_e) = %.3e\nNumber of iterations: %d',xf(6),iter);
    [te, ~, xx_F] = propagate(0,xx0_new, te,mu,false);
    [te, ~, xx_B] = propagate(0,xx0_new, -te,mu,false);
    figure(4)
    plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),LineWidth=2)
    hold on
    plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),LineWidth=2)
    plot(1-mu,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
    plot(-mu,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
    xlabel('x [-]','FontSize',14)
    ylabel('y [-]','FontSize',14)
    zlabel('z[-]','FontSize',14)
%     legend(["Forward","Backward"],'Location','best',FontSize=14)
end

%% Functions

function [dxdt] = xyzCR3BP_noSTM(~,xx, mu)
    
    % Extract variables
    x  = xx(1);
    y  = xx(2);
    z  = xx(3);
    vx = xx(4);
    vy = xx(5);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2+z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2+z^2);

    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
    dUdz = - (1-mu)/r1^3*z -mu/r2^3*z;

    % Assemble right-hand side
    dxdt = zeros(4,1);

    dxdt(1:3) = xx(4:6);
    dxdt(4)   = dUdx + 2*vy;
    dxdt(5)   = dUdy - 2*vx;
    dxdt(6)   = dUdz;
end

function [tf, xf, xx]  = propagate(t0,x0,tf,mu,varargin)
    if nargin>4
        evtFlag=varargin{1};
    else
        evtFlag=1;
    end

    tof = tf - t0;


    % Append to initial conditions the conditions for the STM
    x0Phi0 = x0;
    
    % Perform integration
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@(x,y) x_axis_crossing(x,y,evtFlag));
    [tt, xx] = ode78(@(t,x) xyzCR3BP_noSTM(t,x,mu), [0 tof], x0Phi0, options_STM);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    tf = tt(end);

end

%Event function: stops the integration whenever it is respected the
%condition imposed: for this case when the x - axis is crossed
function [value, isterminal, direction] = x_axis_crossing(~,xx,isTerminal)
    value = xx(2);
    isterminal = isTerminal;
    direction = 0;
end 

