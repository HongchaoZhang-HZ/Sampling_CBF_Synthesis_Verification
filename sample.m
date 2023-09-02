function X = sample(r, num, sampling_flag)
% Given a spherical set, randomly sample points in it.
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
    phi = 2*pi*rand;

    if sampling_flag
%         X(1,i) = r*0.5*(rand+1)*cos(theta)*sin(phi);
%         X(2,i) = r*0.5*(rand+1)*sin(theta)*sin(phi);
%         X(3,i) = r*0.5*(rand+1)*cos(phi);
        X(1,i) = r*rand*cos(theta)*sin(phi);
        X(2,i) = r*rand*sin(theta)*sin(phi);
        X(3,i) = r*rand*cos(phi);
    else
        X(1,i) = r*cos(theta)*sin(phi);
        X(2,i) = r*sin(theta)*sin(phi);
        X(3,i) = r*cos(phi);
    end
end

end