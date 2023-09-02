function [is_invariant,Q] = verify_invariance(A,B,P,invariance_verify_samples,n,D)

cvx_begin sdp

variable Q(n+1,n+1) symmetric
minimize(1)
Q == semidefinite(n+1);
for i=1:D
    -invariance_verify_samples(:,i)'*(P*A + A'*P)*invariance_verify_samples(:,i) == [invariance_verify_samples(:,i) ; 1]'*Q*[invariance_verify_samples(:,i) ; 1];
end

cvx_end

if cvx_optval==1
    is_invariant =1;
else
    is_invariant = 0;
end


