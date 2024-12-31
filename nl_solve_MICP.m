function [Q1_opt, Q2_opt, alpha_opt, b_coeff_opt] = nl_solve_MICP(n, m, d, ...
    M, N, Y_sam, X_sam, Y, f_a_affine, g_a_affine, x, z)
% Solve the MICP to compute the barrier function candidate.
%%%%%%%%%%
% n: the dimension of sample points
% m: the dimension of control inputs
% d: the highest order of the monomials of b(x)
% M: the number of sample points in V[h]
% N: the number of sample points in h(x) >= 0
% Y_sam: a poised set of sample points with h(x) = 0
% X_sam: a poised set of sample points with h(x) >= 0
% Y: a sufficiently large number
% f: the drift term of the system model
% g: the control term of the system model
% x_ex: the external sdp variable of the system state
%%%%%%%%%%
% Q1: the coefficient matrix of the SOS of -b evaluated at V[h]
% Q2: the coefficient matrix of the SOS of dbdx*f evaluated at h(x) >= 0
% alpha_i: binary variables
% b_coeff: the coefficients of the barrier function candidate

% yalmip('clear')

% Define variables
% x = sdpvar(n, 1);
mon_b = monolist(x, d);
mon_sos = monolist(x, d/2);
num_mon_b = size(mon_b, 1);
num_mon_sos = size(mon_sos, 1);
Q1 = sdpvar(num_mon_sos);
Q2 = sdpvar(num_mon_sos);
alpha = binvar(N, 1);
b_coeff = sdpvar(num_mon_b, 1);

% Define objective function
objective = -sum(alpha);

% Define constraints
% constraints = [Q1>=0,Q2>=0];
constraints = [];
for i = 1:num_mon_sos
    constraints = [constraints, Q1(i,i) >= sum(abs(Q1(i,:))) - abs(Q1(i,i)), ...
        Q2(i,i) >= sum(abs(Q2(i,:))) - abs(Q2(i,i))];
end

dbdx = jacobian(b_coeff'*mon_b,x);
dbdxa = dbdx; %db/dz
dbdxa(3) = dbdxa(3)*z(2); %db/dz*dz/dtheta
dbdxg = dbdxa * g_a_affine;
% dbdxg = dbdx * g_a_affine;

% constraints = [constraints, replace(b_coeff'*mon_b, x', [0; 0; 0]') >= 1e-5];
constraints = [constraints, z'*z==1, z<=[1,1]', z>=[-1,-1]', ...
               b_coeff(2)==0, b_coeff(3)==0, b_coeff(4)==0, ...
               b_coeff(5)==0, b_coeff(6)==0, b_coeff(7)==0];
% constraints = [constraints, b_coeff'*mon_b==0, ...
%                dbdx* g_a_affine==[0,0], ...
%                b_coeff(2)==0, b_coeff(3)==0, b_coeff(4)==0, ...
%                b_coeff(5)==0, b_coeff(6)==0, b_coeff(7)==0];


% Set parameters to avoid trivial solutions
epsilon = 1e-3;
for i = 1:M
    constraints = [constraints, replace(b_coeff'*mon_b + mon_sos'*Q1*mon_sos, x', Y_sam(:,i)') + epsilon == 0];
end

for i = 1:N
    constraints = [constraints, replace(b_coeff'*mon_b + Y*abs((1 - alpha(i))), x', X_sam(:,i)') >= 0]; % b==0
    constraints = [constraints, replace(b_coeff'*mon_b - Y*abs((1 - alpha(i))), x', X_sam(:,i)') <= 0]; % b==0
    constraints = [constraints, replace(dbdxg + Y*abs(1 - alpha(i))*ones(1,m), x', X_sam(:,i)') >= 0];
    constraints = [constraints, replace(dbdxg - Y*abs(1 - alpha(i))*ones(1,m), x', X_sam(:,i)') <= 0];
    constraints = [constraints, replace(dbdxa*f_a_affine - mon_sos'*Q2*mon_sos + Y*abs(1 - alpha(i)), x', X_sam(:,i)') >= 0];
    constraints = [constraints, replace(dbdxa*f_a_affine - mon_sos'*Q2*mon_sos - Y*abs(1 - alpha(i)), x', X_sam(:,i)') <= 0];
end

% % Set options for YALMIP solver
options = sdpsettings('solver', 'gurobi', 'gurobi.Threads', 10, 'gurobi.TimeLimit', 120);

% Solve the problem
optimize(constraints, objective, options);

% Retrieve optimal values
Q1_opt = value(Q1);
Q2_opt = value(Q2);
alpha_opt = value(alpha);
% disp(sum(alpha_opt))
b_coeff_opt = value(b_coeff)
f_opt = value(objective);

end