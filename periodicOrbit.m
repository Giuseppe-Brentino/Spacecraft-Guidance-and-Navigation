% Spacecraft Guidance and Navigation (2023/24)
% Assignment 1, exercise 1
% Author: Giuseppe Brentino

clearvars; close all; clc;
plotStyle;

%% Ex 1

mu = 0.012150;
guess = 0;
x_vect = -2:0.001:2;

[L1x,L1y,dUdx] = findLagrangePoint(mu,guess,x_vect);

% plot

figure()
hold on
grid on
plot(x_vect,dUdx,'DisplayName','Gradient of the potential U')
ylim([-30 30])
plot(L1x,L1y,'o','MarkerFaceColor','r','DisplayName','L1 Point')
xlabel('$x$')
ylabel('$\displaystyle\frac{dU}{dx}$')
legend;

%% Ex 2
clearvars; close all; clc;

mu = 0.012150;

% Initial conditions
x0  = 1.08892819445324;
y0  = 0;
z0  = 0.0591799623455459;
vx0 = 0;
vy0 = 0.257888699435051;
vz0 = 0;
phi0 = reshape(eye(6),36,1);

xx0 = [x0;y0;z0;vx0;vy0;vz0;phi0];

%%% Correction of initial states and correct orbit propagation:

% Define solver's parameters
Nmax    = 100;      % Maximum number of iterations
tol     = 1e-14;    % Desired tolerance

% Find correct initial states and propagate the orbit
opt = odeset('RelTol',2.22045e-14,'AbsTol',1e-15,'Events',@(t,xx)planeCrossingEvent(t,xx));
[states, errors] = findOrbit(xx0,tol,Nmax,mu,opt);

%%% plot the required Halo orbit orbit

guess = 2;
[L2x,L2y] = findLagrangePoint(mu,guess);    %compute Lagrange Point L2

figure()
hold on
grid on
plot3(1-mu,0,0,'o','MarkerSize',17,'MarkerFaceColor',[0.5 0.5 0.5],...
    'DisplayName','Moon');
plot3(L2x,L2y,0,'o','MarkerSize',8,'MarkerFaceColor','r',...
    'DisplayName','L2 Point')
plot3(states(:,1),states(:,2),states(:,3),'Color',"#0072BD",...
    'DisplayName','Halo orbit');
plot3(states(:,1),-states(:,2),states(:,3),'Color',"#0072BD",...
    'HandleVisibility','off');
axis equal
view(10,5)
xlabel('x')
ylabel('y')
zlabel('z')
legend

%% Ex 3
clearvars; close all; clc;

mu = 0.012150;
guess = 2;
[L2x,L2y] = findLagrangePoint(mu,guess);    %compute Lagrange Point L2

% definve vector of colors for plot
Color = ["#0072BD"
    "#D95319"
    "#EDB120"
    "#7E2F8E"
    "#77AC30"
    "#4DBEEE"
    "#A2142F"
    "#0080FF"
    "#FFA040"
    "#E57373"
    ];

% Initial conditions
x0  = 1.08892819445324;
y0  = 0;
z0  = 0.0591799623455459;
vx0 = 0;
vy0 = 0.257888699435051;
vz0 = 0;
phi0 = reshape(eye(6),36,1);
xx0 = [x0;y0;z0;vx0;vy0;vz0;phi0];

zz0 = linspace(z0,0.034,10);
z0_length = length(zz0);

% zero-finding algorithm parameters
tol     = 1e-14; % Desired tolerance
Nmax    = 100;   % Maximum number of iterations
opt = odeset('RelTol',2.22045e-14,'AbsTol',1e-15,'Events',@(t,xx)planeCrossingEvent(t,xx));

% start plot
figure()
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
plot3(1-mu,0,0,'o','MarkerSize',17,'MarkerFaceColor',[0.5 0.5 0.5],...
    'DisplayName','Moon');
plot3(L2x,L2y,0,'o','MarkerSize',8,'MarkerFaceColor','r',...
    'DisplayName','L2 Point')

% compute 10 orbits with z0 between the two provided values exploiting
% numerical continuation

for i = 1:z0_length

    if i > 1
        xx0([1;5;3]) = [xx(1,[1 5])';zz0(i)];
    end

    [xx, errors] = findOrbit(xx0,tol,Nmax,mu,opt);
    xx0 = xx(1,:)';

    plot3(xx(:,1),xx(:,2),xx(:,3),'Color',Color(i),'HandleVisibility','off');
    plot3(xx(:,1),-xx(:,2),xx(:,3),'Color',Color(i),'HandleVisibility','off');
    legend;

end
view(10,5)

%% functions

function dx = CRTBP(~,xx,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to integrate the dynamics of a body using the
% Circular Restricted Three Body Problem approximation. It also compute the
% State Transition Matrix.
%
% INPUT:
% - t[1]: time instant. It must be the first input in order to use
%       this function inside an ode solver, but its not required in this
%       case.
% - xx[42,1]: State vector.
% - mu[1]: Gravitational constant of the CRTBP.
%
% OUTPUT:
% -dx[42,1]: time derivative of the state vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s/c states
x = xx(1);
y = xx(2);
z = xx(3);
vx = xx(4);
vy = xx(5);
vz = xx(6);

% Compute distances from bodies 1 and 2
r1 = sqrt((x + mu)^2 + y^2+z^2);
r2 = sqrt((x + mu - 1)^2 + y^2+z^2);
% compute partial derivatives of scalar potential function U
dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
dUdz = - (1-mu)/r1^3*z -mu/r2^3*z;

% assemble right-hand side
dx = zeros(42,1);
dx(1:3) = [vx;vy;vz];
dx(4) = 2*vy + dUdx;
dx(5) = -2*vx + dUdy;
dx(6) = dUdz;


%%% Integration of the STM

% Put the STM PHI in matrix form
phi = reshape(xx(7:end),6,6);

% Jacobian of f
A = zeros(6);
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,1) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1;
A(4,2) = (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(4,3) = (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(4,5) = 2;
A(5,1) = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
A(5,2) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1;
A(5,3) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(5,4) = -2;
A(6,1) = (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
A(6,2) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(6,3) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);

% Compute the derivative of the STM
dphi = A*phi;

dx(7:end) = reshape(dphi,36,1);

end

function [value,isterminal,direction] = planeCrossingEvent(~,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a standard Event function to be used in an ode integrator to stop
% the integration when the second state reaches a value of zero. In this
% framework it is used to stop integrating when the y-plane is crossed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
value = x(2);
isterminal = 1;
direction = 0;
end

function [Lx,Ly,varargout] = findLagrangePoint (mu,guess,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This functions is used to find the three collinear Lagrange Libration
%Points in the 2DCRTBP.
%
% INPUT:
% -mu[1]: Gravitational constant of the CRTBP.
% -guess[1]: Initial guess of the x coordinate of the Lagrange Point.
% -varargin[1]: Vector of x used to evaluate the derivative of the
%               Potential function.
%
% OUTPUT:
%
% -Lx[1]: x-coordinate of the Lagrangian point found.
% -Ly[2]: y-coordinate of the Lagrangian point found.
% -varargout[1]: Derivative of the potential function, evaluated at
%                varargin{1}.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dUdx = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x+mu-1)./abs(x+mu-1).^3 ;
opt = optimoptions('fsolve','OptimalityTolerance',1e-10,'Display','iter-detailed'); % guarantees 10-digit accuracy
Lx = fsolve(dUdx,guess,opt);    %compute Lagrange Point L1
Ly = dUdx(Lx);

if nargout == 3 && nargin == 3
    varargout{1} = dUdx(varargin{1});
end

end

function [states, errors] = findOrbit(xx0,tol,Nmax,mu,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements Newton's method to find the correct initial
% state that guarantees the periodicity of the orbit.
%
% INPUT:
% -xx0[41,1]: Guess of the initial state vector.
% -tol[1]: Tolerance on the error on the initial v_x and v_z.
% -Nmax[1]: Maximum number of iterations.
% -mu[1]:  Gravitational constant of the CRTBP.
% -opt[1]: Structure containing the options of the ode integrator.
%
% OUTPUT:
% -states[n,42]: n-by-42 matrix containing the propagated orbital position
%                velocity and the realtive State Transition Matrix. n is a
%                variable number.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


iter = 0;           % Iteration counter
errors = ones(2,1); % Initial guess of error on vx an vz (higher than tolerance)

% start Newton Method
while (abs(errors(1))>tol || abs(errors(2))>tol) && iter < Nmax

    %%% Compute the new initial conditions

    if iter % if this is not the first iteration, compute the new guesses

        % Select the only elements of the STM useful to correct the previous
        % guess
        phi = reshape(final_states(7:end),6,6);
        phi_p = [phi(4,1), phi(4,5);phi(6,1), phi(6,5)];

        % Compute the new guess
        guess = xx0([1;5]) - phi_p\[final_states(4);final_states(6)];
        xx0([1;5]) = guess;

    end

    % Propagate the orbit with the new initial conditions
    [~,states] = ode113(@(t,xx)CRTBP(t,xx,mu),[0,100],xx0,opt);
    final_states = states(end,:);

    % Compute the error between the reference vx_final and vz_final (=0)
    % and the ones obtained propagating the dynamics
    errors = [final_states(4); final_states(6)];

    %increase iteration counter
    iter = iter + 1;

end

if iter > Nmax
    error("The algorithm didn't converge")
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
set(0,'defaultLineMarkerEdgeColor','k')
set(0,'defaultLineMarkerFaceColor','auto')
% legend:
set(0, 'defaultLegendLocation','best');
set(0, 'defaultLegendFontSize',12);
% axes:
set(0,'defaultAxesFontSize',16);
end
