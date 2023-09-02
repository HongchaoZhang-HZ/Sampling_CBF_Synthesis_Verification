% Nonlinear Verification
close all; clc; clear;
% Define Quadrotor drone's system dynamics
tic
% define constants
D = 128;
% D = 2;
% D = 16; 
r = 1;
[M, N, Y_sam] = initialize(D, r);

% N = 32;
% N = 64;
n = 3;
m = 1;
% d = 4;
d = 2;
Y = 1e8;

x_ex = sdpvar(n, 1);
tau = 0.2;
A = [0, 1, 0;
     0, 0, 1;
     0, 0, -1/tau];
B = [0; 0; 1/tau];
f = A*x_ex;
g = B;

i = 1;

X_sam = sample(r, N, 1);

alpha_opt = zeros(N, 1);

figure(1)
scatter3(X_sam(1,:), X_sam(2,:), X_sam(3,:))
hold on

evalues_CBF_check_st = [];

while sum(alpha_opt) < D_invariance
% is_contained = 0;
% is_invariant = 0;
% while is_contained==0||is_invariant==0
    disp(i)
    
    N = size(X_sam, 2);
    [Q1_opt, Q2_opt, alpha_opt, b_coeff_opt] = nl_solve_MICP(n, m, d, M, ...
        N, Y_sam, X_sam, Y, f_a_affine, g_a_affine, x, z);
%     [Q1_opt, Q2_opt, alpha_opt, b_coeff_opt] = solve_MICP(n, m, d, M, ...
%         N, Y_sam, X_sam, Y, f_l_affine, g_l_affine);
    [beta, P] = decompose_b(b_coeff_opt);
    [is_contained,is_invariant] = verifier_decomposed_bx(beta,P);
    if is_contained && is_invariant
        break;
    else
        X_sam = sample6d([1, 1, 0.1, 1, 1, 1], D+1);
    end
    
    i = i + 1;
end
toc
% Q1_opt
% Q2_opt
% alpha_opt
% b_coeff_opt
syms x1 x2 x3 x4 x5 x6
% b_coeff_opt'*monolist([x1, x2], d)
b_coeff_opt'*monolist([x1 x2 x3 x4 x5 x6], d)
[beta, P] = decompose_b(b_coeff_opt)
[flag_con,flag_inv]=verifier_decomposed_bx(beta, P)