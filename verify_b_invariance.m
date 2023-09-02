function sol = verify_b_invariance(b_coeff, degree, dd, A, g)
% Given a barrier function, verify if it is valid.
%%%%%%%%%%
% b_coeff: the coefficient vector of b(x)
% degree: the order of b(x)
% f: the drift term of the system model
% g: the control term of the system model
%%%%%%%%%%
% X: the set of the sample points


% x = sdpvar(n, 1);
% mon_b = monolist(x, d);
% % syms x1 x2 x3
% % x = [x1 x2 x3]
% b = b_coeff'*mon_b;
% 
% dbdx = diff(b, x);
% 
% mon_eta_theta = monolist(x, dd);
% num_mon_e_t = size(mon_eta_theta, 1);
% eta_coeff = sdpvar(num_mon_e_t, 1);

x1 = sdpvar(1, 1);
x2 = sdpvar(1, 1);
x3 = sdpvar(1, 1);

x = [x1; x2; x3];

% size(A)
% size(x)
f = A*x;

mon_b = monolist(x', degree);
mon_e_t = monolist(x', dd);
% size(mon_b)
% size(b_coeff)
% mon_v = monolist(x', degree/2);
c_eta = sdpvar(length(mon_e_t), 1);
eta = c_eta'*mon_e_t;
c_theta = sdpvar(length(mon_e_t), 1);
theta = c_theta'*mon_e_t;
c_alpha = sdpvar(length(mon_e_t), 1);
alpha = c_alpha'*mon_e_t;

b = b_coeff'*mon_b;
% b = 1122422696836391/9007199254740992 - (1174273219508785*x2)/9903520314283042199192993792 - (729043727400323*x3)/2535301200456458802993406410752 - (4929021559008943*x1*x2)/1125899906842624 - (2292172359870331*x1*x3)/2251799813685248 - (8251820495533165*x2*x3)/4503599627370496 - (1635734464153497*x1^2)/562949953421312 - (8673036083286799*x2^2)/2251799813685248 - (372478008478927*x3^2)/562949953421312 - (6416188664018213*x1)/316912650057057350374175801344;
% beta = 0.086866365688659;
% P = [1.81666666666667	1.15000000000000	0.0833333333333336;
% 1.15000000000000	2.00833333333333	0.150000000000000;
% 0.0833333333333336	0.150000000000000	0.108333333333333];
% b = beta - [x1; x2; x3]'*P*[x1; x2; x3];
% b = x2*(x1/2 - (6753468558286395*x2)/9007199254740992 + x3/2) + x1*(x2/2 - (4971228397594833*x1)/2251799813685248 + (3376734279143197*x3)/4503599627370496) + x3*((3376734279143197*x1)/4503599627370496 + x2/2 - (2533590600238589*x3)/2251799813685248) + 1/20;
% b = x1*x2 + (1688367139571597*x1*x3)/1125899906842624 + x2*x3 - (1805757052820021*x1^2)/562949953421312 - (3940166953256851*x2^2)/2251799813685248 - (4785390413923841*x3^2)/2251799813685248 + 1;
% b = b_coeff'*mon_b;

p = -alpha*jacobian(b, x)*f + eta*b + theta*jacobian(b, x)*g + 1;
% p = jacobian(b, x)*f - eta*b - theta*jacobian(b, x)*g - 1;
% p = -p;
F = [sos(p), sos(alpha)];
% options.solver = 'gurobi';
options = sdpsettings('solver', 'sdpt3', 'verbose', 2);
% options = sdpsettings('solver', 'mosek', 'verbose', 2);
sol = solvesos(F, [c_eta, c_theta, c_alpha], options);
% sol = solvesos(F, [c_eta, c_theta], options);

% % x = sdpvar(n, 1);
% % pvar x1 x2 x3;
% x = [x1; x2; x3];
% 
% prog = sosprogram(x);
% 
% [prog, eta] = sospolyvar(prog,monomials(x,1:degree),'wscoeff');
% [prog, theta] = sospolyvar(prog,monomials(x,1:degree),'wscoeff');
% 
% mon_b = monolist(x, degree);
% b = b_coeff'*mon_b;
% 
% tau = 0.2;
% A = [0, 1, 0;
%      0, 0, 1;
%      0, 0, -1/tau];
% f = A*x;
% dbdxf = diff(b, x)*f;
% dbdxg = diff(b, x)*g;
% 
% % prog = sosineq(prog, dbdxf + eta*b + theta*dbdxg);
% 
% 
% 
% options.solver = 'mosek';
% 
% % [prog, info] = sossolve(prog, options);
% % UNSAT = info.pinf;
% 
% 
% p = dbdxf + eta*b + theta*dbdxg;
% [Q, Z] = findsos(p, options);
end