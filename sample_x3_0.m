function X = sample_x3_0(r, num)
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

if num == 1
    theta = 2*pi*rand;

    X(1,1) = r*0.5*(rand+1)*cos(theta);
    X(2,1) = r*0.5*(rand+1)*sin(theta);

    return
end

for i = 1:num/4
    theta = 2*pi*rand;
    phi = 2*pi*rand;
    
    X(1,i) = r*0.5*(rand+1)*cos(theta)*sin(phi);
    X(2,i) = r*0.5*(rand+1)*sin(theta)*sin(phi);
    X(3,i) = r*0.5*(rand+1)*cos(phi);
    
end

for i = (num/4 + 1):num
    theta = 2*pi*rand;

    X(1,i) = r*0.5*(rand+1)*cos(theta);
    X(2,i) = r*0.5*(rand+1)*sin(theta);

end



end