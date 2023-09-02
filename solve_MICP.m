function [Q1_opt, Q2_opt, alpha_opt, b_coeff_opt] = solve_MICP(n, m, d, M, N, Y_sam, X_sam, Y, f, g, x_ex)
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
% num_mon = 0;
% for i = 1:d
%     num_mon = num_mon + size(combntns(1:n, i), 1);
% end

% x = sdpvar(n-1, 1);
x = sdpvar(n, 1);
mon_b = monolist(x, d);
mon_sos = monolist(x, d/2);
num_mon_b = size(mon_b, 1);
num_mon_sos = size(mon_sos, 1);
Q1 = sdpvar(num_mon_sos);
Q2 = sdpvar(num_mon_sos);
alpha = binvar(N, 1);
b_coeff = sdpvar(num_mon_b, 1);

% Define objective function
% objective = -sum(alpha) - 0.01*norm(b_coeff, 1);
objective = -sum(alpha);
% objective = 0;
% for i = 1:N
%     objective = objective + alpha(i);
% end

% monolist(x, [1, 2])
% b_coeff
% Q1
% 
% constraints = [constraints, append];

% b_coeff'*mon_b+mon_sos'*Q1*mon_sos
% f = replace(mon_b, x, Y(:,1))
% sdisplay(f)
% subs(b_coeff'*mon_b + mon_sos'*Q1*mon_sos, {x(1), x(2), x(3)}, {Y(1,1), Y(2,1), Y(3,1)})

% Define constraints
% constraints = [Q1 >= 0, ...
%                Q2 >= 0, ...
% %                subs(monolist(x, [1, 2])'*b_coeff + monolist(x, [1, 2])'*Q1*monolist(x, [1, 2]) == 0, x, Y)];
%                subs(monolist(x, d)'*b_coeff + monolist(x, d)'*Q1*monolist(x, d) == 0, {x(1), x(2), x(3)}, {[Y(1,1), Y(1,2)], [Y(2,1), Y(2,2)], [Y(3,1), Y(3,2)]})]; 

constraints = [];

% constraints = [Q1 - 0.2*ones(num_mon_sos) >= 0, Q2 - 0.2*ones(num_mon_sos) >= 0];
% constraints = [Q1 - 0.5*eye(num_mon_sos) >= 0, Q2 - 0.5*eye(num_mon_sos) >= 0];

% PSD
% constraints = [Q1 >= 0, Q2 >= 0];

% diagonal dominant
% for i = 1:num_mon_sos
% %     constraints = [constraints, Q1(i,i) >= sum(abs(Q1(i,:))) - abs(Q1(i,i))];
% %     constraints = [constraints, Q2(i,i) >= sum(abs(Q2(i,:))) - abs(Q2(i,i))];
% %     constraints = [constraints, Q1(i,i) >= sum(abs(Q1(i,:))) - abs(Q1(i,i)), Q2(i,i) >= sum(abs(Q2(i,:))) - abs(Q2(i,i))];
% %     constraints = [constraints, Q1(i,i) >= abs(sum(Q1(i,:) - Q1(i,i))), Q2(i,i) >= abs(sum(Q2(i,:) - Q2(i,i)))];
% end

% SOCP
% Q1
M1_12_11 = sdpvar(1, 1);
M1_12_12 = sdpvar(1, 1);
M1_12_21 = sdpvar(1, 1);
M1_12_22 = sdpvar(1, 1);

M1_13_11 = sdpvar(1, 1);
M1_13_13 = sdpvar(1, 1);
M1_13_31 = sdpvar(1, 1);
M1_13_33 = sdpvar(1, 1);

M1_14_11 = sdpvar(1, 1);
M1_14_14 = sdpvar(1, 1);
M1_14_41 = sdpvar(1, 1);
M1_14_44 = sdpvar(1, 1);

M1_23_22 = sdpvar(1, 1);
M1_23_23 = sdpvar(1, 1);
M1_23_32 = sdpvar(1, 1);
M1_23_33 = sdpvar(1, 1);

M1_24_22 = sdpvar(1, 1);
M1_24_24 = sdpvar(1, 1);
M1_24_42 = sdpvar(1, 1);
M1_24_44 = sdpvar(1, 1);

M1_34_33 = sdpvar(1, 1);
M1_34_34 = sdpvar(1, 1);
M1_34_43 = sdpvar(1, 1);
M1_34_44 = sdpvar(1, 1);

constraints = [constraints, Q1(1,1) == M1_12_11 + M1_13_11 + M1_14_11, Q1(1,2) == M1_12_12, Q1(1,3) == M1_13_13, Q1(1,4) == M1_14_14];

constraints = [constraints, Q1(2,1) == M1_12_21, Q1(2,2) == M1_12_22 + M1_23_22 + M1_24_22, Q1(2,3) == M1_23_23, Q1(2,4) == M1_24_24];

constraints = [constraints, Q1(3,1) == M1_13_31, Q1(3,2) == M1_23_32, Q1(3,3) == M1_13_33 + M1_23_33 + M1_34_33, Q1(3,4) == M1_34_34];

constraints = [constraints, Q1(4,1) == M1_14_41, Q1(4,2) == M1_24_42, Q1(4,3) == M1_34_43, Q1(4,4) == M1_14_44 + M1_24_44 + M1_34_44];
% Q1(2,1) = M1_12_21;
% Q1(2,2) = M1_12_22 + M1_23_22 + M1_24_22;
% Q1(2,3) = M1_23_23;
% Q1(2,4) = M1_24_24;
% 
% Q1(3,1) = M1_13_31;
% Q1(3,2) = M1_23_32;
% Q1(3,3) = M1_13_33 + M1_23_33 + M1_34_33;
% Q1(3,4) = M1_34_34;

% Q1(3,1) = M1_13_31;
% Q1(3,2) = M1_23_32;
% Q1(3,3) = M1_13_33 + M1_23_33 + M1_34_33;
% Q1(3,4) = M1_34_34;

% Q1(4,1) = M1_14_41;
% Q1(4,2) = M1_24_42;
% Q1(4,3) = M1_34_43;
% Q1(4,4) = M1_14_44 + M1_24_44 + M1_34_44;

constraints = [constraints, M1_12_11 + M1_12_22 >= 0];
constraints = [constraints, M1_13_11 + M1_13_33 >= 0];
constraints = [constraints, M1_14_11 + M1_14_44 >= 0];
constraints = [constraints, M1_23_22 + M1_23_33 >= 0];
constraints = [constraints, M1_24_22 + M1_24_44 >= 0];
constraints = [constraints, M1_34_33 + M1_34_44 >= 0];

constraints = [constraints, cone([2*M1_12_12; M1_12_11 - M1_12_22], M1_12_11 + M1_12_22)];
constraints = [constraints, cone([2*M1_13_13; M1_13_11 - M1_13_33], M1_13_11 + M1_13_33)];
constraints = [constraints, cone([2*M1_14_14; M1_14_11 - M1_14_44], M1_14_11 + M1_14_44)];
constraints = [constraints, cone([2*M1_23_23; M1_23_22 - M1_23_33], M1_23_22 + M1_23_33)];
constraints = [constraints, cone([2*M1_24_24; M1_24_22 - M1_24_44], M1_24_22 + M1_24_44)];
constraints = [constraints, cone([2*M1_34_34; M1_34_33 - M1_34_44], M1_34_33 + M1_34_44)];
% 
% % Q2
M2_12_11 = sdpvar(1, 1);
M2_12_12 = sdpvar(1, 1);
M2_12_21 = sdpvar(1, 1);
M2_12_22 = sdpvar(1, 1);

M2_13_11 = sdpvar(1, 1);
M2_13_13 = sdpvar(1, 1);
M2_13_31 = sdpvar(1, 1);
M2_13_33 = sdpvar(1, 1);

M2_14_11 = sdpvar(1, 1);
M2_14_14 = sdpvar(1, 1);
M2_14_41 = sdpvar(1, 1);
M2_14_44 = sdpvar(1, 1);

M2_23_22 = sdpvar(1, 1);
M2_23_23 = sdpvar(1, 1);
M2_23_32 = sdpvar(1, 1);
M2_23_33 = sdpvar(1, 1);

M2_24_22 = sdpvar(1, 1);
M2_24_24 = sdpvar(1, 1);
M2_24_42 = sdpvar(1, 1);
M2_24_44 = sdpvar(1, 1);

M2_34_33 = sdpvar(1, 1);
M2_34_34 = sdpvar(1, 1);
M2_34_43 = sdpvar(1, 1);
M2_34_44 = sdpvar(1, 1);

constraints = [constraints, Q2(1,1) == M2_12_11 + M2_13_11 + M2_14_11, Q2(1,2) == M2_12_12, Q2(1,3) == M2_13_13, Q2(1,4) == M2_14_14];

constraints = [constraints, Q2(2,1) == M2_12_21, Q2(2,2) == M2_12_22 + M2_23_22 + M2_24_22, Q2(2,3) == M2_23_23, Q2(2,4) == M2_24_24];

constraints = [constraints, Q2(3,1) == M2_13_31, Q2(3,2) == M2_23_32, Q2(3,3) == M2_13_33 + M2_23_33 + M2_34_33, Q2(3,4) == M2_34_34];

constraints = [constraints, Q2(4,1) == M2_14_41, Q2(4,2) == M2_24_42, Q2(4,3) == M2_34_43, Q2(4,4) == M2_14_44 + M2_24_44 + M2_34_44];

constraints = [constraints, M2_12_11 + M2_12_22 >= 0];
constraints = [constraints, M2_13_11 + M2_13_33 >= 0];
constraints = [constraints, M2_14_11 + M2_14_44 >= 0];
constraints = [constraints, M2_23_22 + M2_23_33 >= 0];
constraints = [constraints, M2_24_22 + M2_24_44 >= 0];
constraints = [constraints, M2_34_33 + M2_34_44 >= 0];

constraints = [constraints, cone([2*M2_12_12; M2_12_11 - M2_12_22], M2_12_11 + M2_12_22)];
constraints = [constraints, cone([2*M2_13_13; M2_13_11 - M2_13_33], M2_13_11 + M2_13_33)];
constraints = [constraints, cone([2*M2_14_14; M2_14_11 - M2_14_44], M2_14_11 + M2_14_44)];
constraints = [constraints, cone([2*M2_23_23; M2_23_22 - M2_23_33], M2_23_22 + M2_23_33)];
constraints = [constraints, cone([2*M2_24_24; M2_24_22 - M2_24_44], M2_24_22 + M2_24_44)];
constraints = [constraints, cone([2*M2_34_34; M2_34_33 - M2_34_44], M2_34_33 + M2_34_44)];


constraints = [constraints, b_coeff(2) == 0, b_coeff(3) == 0, b_coeff(4) == 0];

% eps3 = 1e-9;
eps3 = 0;
% eps3 = 1e-1;
constraints = [constraints, replace(b_coeff'*mon_b, x', [0; 0; 0]') >= eps3];
% constraints = [constraints, b_coeff(1) >= eps3];

epsilon = 1e-2;
% epsilon = 0;
% epsilon = 0.5*1e-3;
eps2 = 1e-5;
for i = 1:M
%     constraints = [constraints, replace(b_coeff'*mon_b + mon_sos'*Q1*mon_sos, x', Y_sam(:,i)') + b_coeff(1) == 0];
%     constraints = [constraints, replace((b_coeff.*[1e3; ones(num_mon_b-1,1)])'*mon_b + mon_sos'*Q1*mon_sos, x', Y_sam(:,i)') + epsilon == 0];
%     constraints = [constraints, replace(b_coeff'*mon_b + mon_sos'*Q1*mon_sos, x', Y_sam(1:2,i)') == 0];
    constraints = [constraints, replace(b_coeff'*mon_b + mon_sos'*Q1*mon_sos, x', Y_sam(:,i)') + epsilon == 0];
end

% jacobian(b_coeff'*mon_b, x_ex)*g 
% Y*(1 - alpha(i))*ones(m,1)
% replace(jacobian(b_coeff'*mon_b, x)*f - mon_sos'*Q2*mon_sos + Y*(1 - alpha(1))*ones(m,1), [x', x_ex'], [X_sam(:,1)', X_sam(:,1)'])

for i = 1:N
%     constraints = [constraints, replace(b_coeff'*mon_b + Y*(1 - alpha(i)), x', X_sam(1:2,i)') >= 0];
%     constraints = [constraints, replace(b_coeff'*mon_b - Y*(1 - alpha(i)), x', X_sam(1:2,i)') <= 0];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x_ex)*g + Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(1:2,i)', X_sam(:,i)']) >= 0];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x_ex)*g - Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(1:2,i)', X_sam(:,i)']) <= 0];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x_ex)*f - mon_sos'*Q2*mon_sos + Y*(1 - alpha(i)) - epsilon, [x', x_ex'], [X_sam(1:2,i)', X_sam(:,i)']) >= 0];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x_ex)*f - mon_sos'*Q2*mon_sos - Y*(1 - alpha(i)) - epsilon, [x', x_ex'], [X_sam(1:2,i)', X_sam(:,i)']) <= 0];
%     replace(b_coeff'*mon_b - Y*abs((1 - alpha(i))), x', X_sam(:,i)')
    constraints = [constraints, replace(b_coeff'*mon_b + Y*abs((1 - alpha(i))), x', X_sam(:,i)') >= 0];
    constraints = [constraints, replace(b_coeff'*mon_b - Y*abs((1 - alpha(i))), x', X_sam(:,i)') <= 0];
    constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*g + Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) >= 0];
    constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*g - Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) <= 0];
    constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*f - mon_sos'*Q2*mon_sos + Y*(1 - alpha(i)), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) >= 0];
    constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*f - mon_sos'*Q2*mon_sos - Y*(1 - alpha(i)), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) <= 0];
%     constraints = [constraints, replace(b_coeff'*mon_b + Y*(1 - alpha(i)), x', X_sam(:,i)') >= 0 - eps2];
%     constraints = [constraints, replace(b_coeff'*mon_b - Y*(1 - alpha(i)), x', X_sam(:,i)') <= 0 + eps2];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*g + Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) >= 0 - eps2];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*g - Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) <= 0 + eps2];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*f - mon_sos'*Q2*mon_sos + Y*(1 - alpha(i)) - epsilon, [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) >= 0];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*f - mon_sos'*Q2*mon_sos - Y*(1 - alpha(i)) - epsilon, [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) <= 0];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*f - mon_sos'*Q2*mon_sos + Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) >= 0];
%     constraints = [constraints, replace(jacobian(b_coeff'*mon_b, x)*f - mon_sos'*Q2*mon_sos - Y*(1 - alpha(i))*ones(m,1), [x', x_ex'], [X_sam(:,i)', X_sam(:,i)']) <= 0];
end

% check(constraints)

% constraints

% Set options for YALMIP solver
options = sdpsettings('solver', 'gurobi');
% options = sdpsettings('solver', 'bnb');
% options = sdpsettings('solver', 'mosek', 'mosek.MSK_IPAR_OPTIMIZER', 'MSK_OPTIMIZER_MIXED_INT_CONIC');
% options = sdpsettings('solver', 'sdpt3');
% options = sdpsettings('solver', 'mosek');
% options.mosek.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_MIXED_INT_CONIC';
% options.mosek.MSK_IPAR_INTPNT_SOLVER = 'MSK_OPTIMIZER_INTPNT';

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