function [M, N, Y] = initialize(D, r)
% Initialize the parameters and samples that will be used in the iterative procedure.
%%%%%%%%%%
% D: 
% r: the radius of the safe region h(x) >= 0
%%%%%%%%%%
% M: the number of sample points in V[h]
% N: the number of sample points in h(x) >= 0
% Y: a poised set of sample points with h(x) = 0

% M = 8*D;
% % M = D/2;
% N = D/2;
% % Y = sample_2d(r, M, 0);
% Y = sample(r, M, 0);

n = 3;

containment_poised = 0;
boundary_samples = [];
D_containment = 0;
while containment_poised==0
    %The sample y should be chosen on the boundary of the safe region. The
    %below code assumes that the boundary points are equal to y'*y = 1.
    y = randn(n,1);
    y_norm = sqrt(y'*y);
    boundary_samples = [boundary_samples y/y_norm];
    containment_poised = check_poisedness(boundary_samples, n+1);
    D_containment = D_containment + 1;
end

% M = 8*D_containment;
M = D_containment;
N = D/2;

% X = zeros(3, 1);
% for i = 1:7*D_containment
%     theta = 2*pi*rand;
%     phi = 2*pi*rand;
% 
%     X(1,1) = r*cos(theta)*sin(phi);
%     X(2,1) = r*sin(theta)*sin(phi);
%     X(3,1) = r*cos(phi);
% 
%     boundary_samples = [boundary_samples X];
% end
Y = boundary_samples;

end