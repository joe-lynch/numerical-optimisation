function [H,H_] = sr1(H,H_,dn,yn,eta)

%     This was tested by comparing it the exact Hessian and exact Hessian
%     inverse for the brockrosen function.
% From tr_dogleg(@f_Rosen,@g_Rosen,eye(2),inv(eye(2)),x0,0.2,1,0.125,1.0e-5)
%     x0 = [-1.2;1]
%     - dogleg + exact hessian : 19 iterations.
%     - dogleg + sr1 update    : 47 iterations.
%     - cauchy + exact hessian : 1209 iterations. 
%     - cauchy + sr1 update    : 1183 iterations.
%     x0 = [1;-1]
%     - dogleg + exact hessian : 10 iterations.
%     - dogleg + sr1 update    : 26 iterations.
%     - cauchy + exact hessian : 807 iterations.
%     - cauchy + sr1 update    : 1246 iterations.
%
% Moreover, calling
%
% [H,H_] = sr1(h_Rosen([-1.2;1]), inv(h_Rosen([-1.2;1])), [0.189;0.0652], [23.246;9.665], 10^(-6))
%
% gives the sr1 update of the hessian for [-1.0109;1.0652]. This can be
% checked with the actual hessian and inverse to see how accurate the update is.



  % negation of update procedure from assignment to avoid else statement
  if not((norm(dn) == 0) || (norm(dn'*(yn - H*dn)) < eta*norm(dn)*norm(yn - H*dn)))
    % compute hessian and then inverse
    H = H + ((yn - H*dn)*(yn - H*dn)')/(dn'*(yn - H*dn));
    H_ = H_ + ((dn - H_*yn)*(dn - H_*yn)')/((dn - H_*yn)'*yn);
  end
end