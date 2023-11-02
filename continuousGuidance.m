% Spacecraft Guidance and Navigation (2023/24)
% Assignment 1, exercise 3
% Author: Giuseppe Brentino

%% Ex 2 
clearvars; close all; clc;
rng default;
addpath('.\kernels')
addpath('.\mice\src\mice')
addpath('.\mice\lib')

cspice_furnsh('kernels\naif0012.tls'); % (LSK)
cspice_furnsh('kernels\de432s.bsp');   % (SPK)
cspice_furnsh('kernels\gm_de432.tpc');  
cspice_furnsh('kernels\pck00010.tpc'); 

t_launch = cspice_str2et('May 28 14:13:09.000 UTC 2023');

m0 = 1000; %kg
Tmax = 800e-6; % KN kg m/s^2  
Isp = 3120; %s
g0 = 9.81e-3; % km/s^2

mu_s = cspice_bodvrd('SUN','GM',1);
l = cspice_convrt(1,'AU','km');
m = m0;
t = sqrt(l^3/(mu_s));

data.center = 'SUN';
data.frame = 'ECLIPJ2000'; 

xi = cspice_spkezr('Earth',t_launch,data.frame,'NONE',data.center);

x0 = [xi;m0];
x0 = [x0(1:3)/l;x0(4:6)*(t/l);x0(7)/m];
ri = x0(1:3);
vi = x0(4:6);

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

opt = optimoptions("fsolve",'Display','iter-detailed');
exitflag = 0;
i = 0;
N_iter = 10;
 fval_vect = zeros(N_iter,1);
 guess_mat = zeros(8,N_iter);

while exitflag ~=1 && i<N_iter
     guess = [40*rand(6,1) - 20;20*rand(1); 2*pi*rand(1)+data.t0];
    [sol,fval,exitflag] = fsolve(@optimization,guess,opt,x0,data);

    i = i+1;
    fval_vect(i) = norm(fval);
    guess_mat(:,i) = sol;
end
fval_vect = fval_vect(fval_vect~=0);
[~,ind] = min(fval_vect);
guess_mat = guess_mat(:,fval_vect~=0);
sol = guess_mat(:,ind);
[T,xx] = ode113(@TPBVP,[data.t0,sol(8)],[x0;sol(1:7)],data.ode_opt,data);
venus = cspice_spkezr('Venus',sol(8)*t,data.frame,'NONE',data.center);

if exitflag == 1
    fprintf('CONVERGEEEEEEEEEE! \n')
    S(1) = load('gong');
    sound(S(1).y,S(1).Fs)
else
    fprintf('RITENTA, sarai più fortunato \n')
    S(1) = load('gong');
S(2) = load('handel');
S(3).y = ((S(2).y(1:length(S(1).y),1)).*(S(1).y));
S(3).Fs = 10000;
sound(S(3).y,S(3).Fs)
 end

 % errori
 err_pos1 = norm(venus(1:3)-xx(end,1:3)'*l);
 err_vel1 = norm(venus(4:6)-xx(end,4:6)'*l/t)*1e3;
 % plot traiettoria
  figure()
   hold on
   grid on
   plot3(xx(:,1),xx(:,2),xx(:,3),'b-','DisplayName','Trajectory')
    plot3(venus(1)/data.l,venus(2)/data.l,venus(3)/data.l,'o','MarkerFaceColor',...
        'r','MarkerSize',5,'DisplayName','Venus at tf')
    plot3(xx(end,1),xx(end,2),xx(end,3),'o','MarkerFaceColor',...
        'c','MarkerSize',5,'DisplayName','S/C at tf')
    legend('Location','best');
    axis equal

    % plot thrust direction
    figure()
    hold on
    grid on
    plot3(xx(:,11),xx(:,12),xx(:,13))
  
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

%% Ex 4 
close all; clc;

opt = optimoptions("fsolve",'Display','iter-detailed');
Tmax = linspace(700,500,5)*1e-6;

for j = 1:length(Tmax)
exitflag = 0;
i = 0;
N_iter = 10;
 fval_vect = zeros(N_iter,1);
 guess_mat = zeros(8,N_iter);
 data.Tmax = Tmax(j)*((t^2)/(m*l));
while exitflag ~=1 && i<N_iter
    guess = sol;
    [sol,fval,exitflag] = fsolve(@optimization,guess,opt,x0,data);
     
    i = i+1;
    fval_vect(i) = norm(fval);
    guess_mat(:,i) = sol;
end
end
fval_vect = fval_vect(fval_vect~=0);
[~,ind] = min(fval_vect);
guess_mat = guess_mat(:,fval_vect~=0);
sol = guess_mat(:,ind);
[T,xx] = ode113(@TPBVP,[data.t0,sol(8)],[x0;sol(1:7)],data.ode_opt,data);
venus = cspice_spkezr('Venus',sol(8)*t,data.frame,'NONE',data.center);

if exitflag == 1
    fprintf('CONVERGEEEEEEEEEE! \n')
    S(1) = load('gong');
    sound(S(1).y,S(1).Fs)
else
    fprintf('RITENTA, sarai più fortunato \n')
 end

 % errori
 err_pos2 = norm(venus(1:3)-xx(end,1:3)'*l);
 err_vel2 = norm(venus(4:6)-xx(end,4:6)'*l/t)*1e3;
 % plot traiettoria
  figure()
   hold on
   grid on
   plot3(xx(:,1),xx(:,2),xx(:,3),'b-','DisplayName','Trajectory')
    plot3(venus(1)/data.l,venus(2)/data.l,venus(3)/data.l,'o','MarkerFaceColor',...
        'r','MarkerSize',5,'DisplayName','Venus at tf')
    plot3(xx(end,1),xx(end,2),xx(end,3),'o','MarkerFaceColor',...
        'c','MarkerSize',5,'DisplayName','S/C at tf')
    legend('Location','best');
    axis equal

    % plot thrust direction
    figure()
    hold on
    grid on
    plot3(xx(:,11),xx(:,12),xx(:,13))

%% functions

function dx = TPBVP(~,x,data)

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
   dyn       = [vv_f;
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