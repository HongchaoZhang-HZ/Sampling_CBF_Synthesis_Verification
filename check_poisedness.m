function poised = check_poisedness(samples, n)

%The below code is based on the test case of a three-dimensional system with quadratic safety constraint and barrier function 
%For this case, the space of polynomials is spanned by u(x)u(x)', where u(x) is the vector of monomials [x1 x2 x3 1]'

% n=4;

[~,num_samples] = size(samples);

U2 = zeros(n^2, num_samples);

for i=1:num_samples
    z = samples(:,i);
    uz = [z ; 1];
    U2(:,i) = vec(uz*uz');
end

if rank(U2) == num_samples
    poised = 0;
else
    poised=1;
end