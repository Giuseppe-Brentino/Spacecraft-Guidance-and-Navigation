% Spacecraft Guidance and Navigation (2023/24)
% Assignment 1, exercise 1
% Author: Giuseppe Brentino

clearvars; close all; clc;

%% Ex 1
clearvars; close all; clc;
mu = 0.012150;
dUdx = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x+mu-1)./abs(x+mu-1).^3 ;
opt = optimoptions('fsolve','OptimalityTolerance',1e-10);
L1 = fsolve(dUdx,0,opt);

x_vect = -2:0.001:2;
figure
hold on
grid on
plot(x_vect,dUdx(x_vect))
ylim([-30 30])

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

% correction of initial states with Pseudo-Newton method
% Initialization
err_vxf = 1;    % Initial guess of error (higher than tolerance)
err_vzf = 1;    % Initial guess of error (higher than tolerance)
Nmax    = 50;   % Maximum number of iterations 
iter    = 0;    % Iteration counter
tol     = 1e-8; % Desired tolerance

% Set the update to the first guess 
vy0_new = vy0;
x0_new  = x0;
vxf_ref = 0;
vzf_ref = 0;
tic
while (abs(err_vxf)>tol || abs(err_vzf)>tol) && iter < Nmax
    % Vector of initial conditions
    xx0 = [x0_new;y0;z0;vx0;vy0_new;vz0;phi0];
    opt = odeset('RelTol',1e-12,'AbsTol',1e-13,'Events',@(t,xx)planeCrossingEvent(t,xx));
    [~,xx] = ode113(@(t,xx)CRTBP(t,xx,mu),[0,100],xx0,opt);
    xx = xx(end,:);

    phi = reshape(xx(7:end),6,6);
    phi_p = [phi(4,1), phi(4,5);phi(6,1), phi(6,5)];

    guess = [x0_new;vy0_new] - phi_p\[xx(4);xx(6)];

    x0_new = guess(1);
    vy0_new = guess(2);
    err_vxf = xx(4);
    err_vzf = xx(6);
    iter = iter + 1;
end

 xx0 = [x0_new;y0;z0;vx0;vy0_new;vz0;phi0];
    opt = odeset('RelTol',1e-12,'AbsTol',1e-13);
    [~,xx] = ode113(@(t,xx)CRTBP(t,xx,mu),[0,15],xx0,opt);

    figure()
    hold on 
    grid on
    plot3(xx(:,1),xx(:,2),xx(:,3));
    plot3(xx(1,1),xx(1,2),xx(1,3),'o','MarkerSize',5,'MarkerFaceColor','r');
    axis equal
    xlabel('x')
    ylabel('y')
%% functions

function dx = CRTBP(~,xx,mu)

% s/c states
dx = zeros(42,1);
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
dx(1:3) = [vx;vy;vz];
dx(4) = 2*vy + dUdx;
dx(5) = -2*vx + dUdy;
dx(6) = dUdz;


%STM
% Put PHI in matrix form
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
    value = x(2);
    isterminal = 1;
    direction = 0;
end













