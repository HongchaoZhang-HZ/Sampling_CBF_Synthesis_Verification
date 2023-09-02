function X = sample_on_sphere(r, num)
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

% if num == 1
%     theta = 2*pi*rand;
% 
%     X(1,1) = 0.9*r*cos(theta);
%     X(2,1) = 0.9*r*sin(theta);
% 
%     return
% end
% 
% 
% for i = 1:num/4
%     theta = 2*pi*rand;
%     phi = 2*pi*rand;
%     
%     X(1,i) = 0.9*r*cos(theta)*sin(phi);
%     X(2,i) = 0.9*r*sin(theta)*sin(phi);
%     X(3,i) = 0.9*r*cos(phi);
%     
% end
% 
% 
% for i = (num/4 + 1):num
%     theta = 2*pi*rand;
% 
%     X(1,i) = 0.9*r*cos(theta);
%     X(2,i) = 0.9*r*sin(theta);
% 
% end

for i = 1:num
    theta = 2*pi*rand;
    phi = 2*pi*rand;

    X(1,i) = r*0.9*cos(theta)*sin(phi);
    X(2,i) = r*0.9*sin(theta)*sin(phi);
    X(3,i) = r*0.9*cos(phi);
end

end