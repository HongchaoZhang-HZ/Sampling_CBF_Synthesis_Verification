% sos synth
% clc; close all; clear all;

pvar x1 x2 x3
x = [x1;x2;x3];
degree = 3;
prog = sosprogram(x);
[prog,beta] = sospolyvar(prog,monomials(x,0:degree),'wscoeff');
% h = -x1^2-x2^2-x3^2+1;

tau = 0.2;
A = [0, 1, 0;
     0, 0, 1;
     0, 0, -1/tau];
B = [0; 0; 1/tau];
g = B;
C = eye(3);
sys = ss(A,B,C,0);
Q = eye(3);
R = 1;
[K,~,~] = lqr(sys,Q,R);
F = A-B*K;
N = eye(3);
P = lyap(F',N);
% b = beta - x'*P*x;
mon_b = monomials(x, 0:2);
b = b_coeff_opt'*mon_b;
% size(A)
% size(x)
f = A*x;

dbdx = jacobian(b, x);
% p = jacobian(b, x)*f + eta*b + theta*jacobian(b, x)*g;

dbdxf = dbdx*f;
dbdxg = dbdx*g;

prog = sosineq(prog, dbdxf);
prog = soseq(prog,dbdxg);
prog = sosineq(prog,b);
prog = sosineq(prog,-b);
% prog = sosineq(prog,h);

solver_opt.solver = 'sdpt3';
[prog,info] = sossolve(prog,solver_opt);
SOLV = sosgetsol(prog,b)