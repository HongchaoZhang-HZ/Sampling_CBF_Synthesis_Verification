%% verify_b with increasing order
dd = 2;
tau = 0.2;
A = [0, 1, 0;
     0, 0, 1;
     0, 0, -1/tau];
g = [0; 0; 1/tau];

% x2*(x1/2 - (1688367139571603*x2)/2251799813685248 + x3/2) + x3*((1688367139571597*x1)/2251799813685248 + x2/2 - (2533590600238593*x3)/2251799813685248) + x1*(x2/2 - (1242807099398709*x1)/562949953421312 + (1688367139571597*x3)/2251799813685248) + 1/100
% b_coeff_opt = [0.01; 0; 0; 0; -2.2077; 1; -0.7498; 1.4996; 1; -1.1251];

while 1
    disp('The current order is:')
    disp(dd)
    sol = verify_b(b_coeff_opt, 2, dd, A, g);
    dd = dd+1;
    if ~(strcmp(sol.info, 'Infeasible problem (MOSEK)')) && ~(strcmp(sol.info, 'Unbounded objective function (MOSEK)'))
        break;
    end
end