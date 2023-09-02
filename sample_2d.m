function X = sample_2d(r, num, sampling_flag)
% Given a circular set, randomly sample points in it.
%%%%%%%%%%
% r: the radius of the spherical set
% num: the number of the sample points
% sampling_flag: whether we sample on the boundary only or in the interior 
%   0: sample on the boundary only
%   1: sample in the interior and the boundary
%%%%%%%%%%
% X: the set of the sample points

n = 3;
X = zeros(n, num);

for i = 1:num
    theta = 2*pi*rand;

    if sampling_flag
        X(1,i) = r*rand*cos(theta);
        X(2,i) = r*rand*sin(theta);
    else
        X(1,i) = r*cos(theta);
        X(2,i) = r*sin(theta);
    end
end

end