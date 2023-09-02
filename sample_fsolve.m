function X_sam = sample_fsolve(N, beta, P, dbdxg)

X_sam = zeros(6, N);
i = 1;
while i <= N % generate 10 random solutions
%     x0 = [7, 8.8882];
%     x0 = 16*(rand(1,2)-0.5); % generate two random numbers between -10 and 10
    x0 = 2*(rand(6, 1) - 0.5);
    [sol, fval, exitflag] = fsolve(@(x) dbdxg{1}, x0);
%     [sol, fval, exitflag] = fsolve(@(x) beta - x'*P*x, x0);
%     [sol, fval, exitflag] = fsolve(@(x)myfun(x, beta, P, dbdxg), x0);
    if exitflag > 0 % if a solution is found
%         disp(sol) % display the solution
%         disp(fval)
%         disp(norm(sol))
        X_sam(:, i) = sol;
        i = i + 1;
    end
%     i = i+1;
end

function dbdxg = myfun(x, beta, P, dbdxg)
% function [b, dbdxg] = myfun(x, beta, P, dbdxg)
    q(1) = x(1);
    q(2) = x(2);
    q(3) = x(3);
    qdot(1) = x(4);
    qdot(2) = x(5);
    qdot(3) = x(6);
%     state = [q(1); q(2); q(3); qdot(1); qdot(2); qdot(3)];
    b = beta - x'*P*x;
%     dbdxg = jacobian(b, x)*g;
%     x1 = x(1);
%     x2 = x(2);
%     x3 = x(3);
%     b = x1*x2 + (1688367139571597*x1*x3)/1125899906842624 + x2*x3 - (1805757052820021*x1^2)/562949953421312 - (3940166953256851*x2^2)/2251799813685248 - (4785390413923841*x3^2)/2251799813685248 + 1;
%     x1 = sdpvar(1);
%     x2 = sdpvar(1);
%     x3 = sdpvar(1);
%     x = [x1, x2, x3];
%     dbdxg = jacobian(b, x);
%     dbdxg = diff(b);
end

% function [b, dbdxg] = myfun(x)
%     x1 = x(1);
%     x2 = x(2);
%     x3 = x(3);
%     b = x1*x2 + (1688367139571597*x1*x3)/1125899906842624 + x2*x3 - (1805757052820021*x1^2)/562949953421312 - (3940166953256851*x2^2)/2251799813685248 - (4785390413923841*x3^2)/2251799813685248 + 1;
%     x1 = sdpvar(1);
%     x2 = sdpvar(1);
%     x3 = sdpvar(1);
%     x = [x1, x2, x3];
%     dbdxg = jacobian(b, x);
% %     dbdxg = diff(b);
% end

end
% function [F, J] = myfun(x)
% F = [x(1)^2 + x(2)^2 - 128;     2*x(1) + 2*x(2)];
% J = [2*x(1), 2*x(2);     2, 2];
% end



% for i = 1:10 % generate 10 random solutions
%     x0 = 20*(rand(1,2)-0.5); % generate two random numbers between -10 and 10
%     [sol, fval, exitflag] = fsolve(@myfun, x0);
%     if exitflag > 0 % if a solution is found
%         disp(sol) % display the solution
%     end
% end
% 
% function [F, J] = myfun(x)
% F = [x(1)^2 + x(2)^2 - 5 + x(1)*x(2);     2*x(1) + x(2);     x(1) + 2*x(2)];
% J = [2*x(1) + x(2), x(1);     2, 1;     1, 2];
% end