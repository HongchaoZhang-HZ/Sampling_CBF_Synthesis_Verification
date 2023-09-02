function X = sample_poisedness(r, sampling_flag)
% Given a spherical set, randomly sample points in it.
%%%%%%%%%%
% r: the radius of the spherical set
% sampling_flag: whether we sample on the boundary only or in the interior 
%   0: sample on the boundary only
%   1: sample in the interior and the boundary
%%%%%%%%%%
% X: the set of the sample points

n = 3;
sample = zeros(n, 1);

invariance_verify_samples = [];
invariance_poised = 0;

while invariance_poised == 0
    theta = 2*pi*rand;
    phi = 2*pi*rand;

    if sampling_flag
        sample(1,1) = r*0.5*(rand+1)*cos(theta)*sin(phi);
        sample(2,1) = r*0.5*(rand+1)*sin(theta)*sin(phi);
        sample(3,1) = r*0.5*(rand+1)*cos(phi);
    else
        sample(1,1) = r*cos(theta)*sin(phi);
        sample(2,1) = r*sin(theta)*sin(phi);
        sample(3,1) = r*cos(phi);
    end
    invariance_verify_samples = [invariance_verify_samples sample];
    invariance_poised = check_poisedness(invariance_verify_samples);
end

X = invariance_verify_samples;

end