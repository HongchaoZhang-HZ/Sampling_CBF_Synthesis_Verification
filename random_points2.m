clc; close all; clear all; %echo on;

%% system model
n = 3;
m = 1;
tau = 0.2;

A = [0, 1, 0;
    0, 0, 1;
    0, 0, -1/tau];

B = [0; 0; 1/tau];

Q = eye(n);
R = eye(m);


%% parameter initialization
T = 1e1;
num_steps = 1e5;
dt = T/num_steps;


%% 1
N = 5;
p = 2;
r_b = 0.5;
x = zeros(n, N);
z = zeros(p+1, N);
z0 = zeros(p+1,1);

for i = 1:N
    theta = 2*pi*rand;
    x(1,i) = r_b*cos(theta);
    x(2,i) = r_b*sin(theta);

%     [sol, fval, exitflag] = fsolve(@(z)[0, 1, -1]*(z.^2), z0);
%     if exitflag > 0 
%         disp(sol) 
%     end
%     [0, 1, -1]*(z.^2) = 0;
end


%% 2
M = 5;
r_h = 1;
y = zeros(n, M);

for i = 1:M
    theta = 2*pi*rand;
    y(1,i) = r_h*cos(theta);
    y(2,i) = r_h*sin(theta);
end


%% 3
pvar x1 x2 x3 z1 z2 z3;
vars = [x1; x2; x3];

f = A*vars;

prog = sosprogram(vars);

x_mon = monomials(vars, 2);
[prog, F_1] = sossosvar(prog, x_mon);
% [prog, F_1] = sospolyvar(prog, x_mon);

b = r_b^2 - x1^2 - x2^2;
dbdx = [diff(b, x1), diff(b, x2), diff(b, x3)];
% dbdx = [-2*x1, -2*x2];

u_max = 1;
u_min = -1;
bu = [u_max; -u_min];

% prog = soseq(prog, dbdx*f);
% prog = soseq(prog, [dbdx*f; bu]);
% prog = soseq(prog, subs([dbdx*f; bu], vars, x));

solver_opt.solver = 'sdpt3';
% [Q, Z] = findsos(F_1, solver_opt);
[prog,info] = sossolve(prog,solver_opt);
% prog = sossolve(prog, solver_opt);

SOLV = sosgetsol(prog, F_1);

% echo off;