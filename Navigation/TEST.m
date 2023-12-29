close all; clearvars; clc;
syms x(t) y(t) z(t) vx(t) vy(t) vz(t) n x0 y0 z0 vx0 vy0 vz0

ode1 = diff(x) == vx;
ode2 = diff(y) == vy;
ode3 = diff(z) == vz;
ode4 = diff(vx) == 3*n^2*x + 2*n*vy;
ode5 = diff(vy) == -2*n*vx;
ode6 = diff(vz) == -n^2*z;

odes = [ode1;ode2;ode3;ode4;ode5;ode6];
cond = [x(0) == x0; y(0) == y0; z(0) == z0; vx(0) == vx0; vy(0) == vy0; vz(0) == vz0;];
Sol = dsolve(odes,cond);
vxd = simplify(rewrite(Sol.vx, "sincos"))
vyd = simplify(rewrite(Sol.vy, "sincos"))
vzd = simplify(rewrite(Sol.vz, "sincos"))