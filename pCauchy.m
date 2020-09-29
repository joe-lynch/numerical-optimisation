function [pC] = pCauchy(H, g, delta)
%PCAUCHY Find the Cauchy Point for the trust region model of f
%   Input:
%       H: Approximated hessian of f at x_k
%       g: Approximated gradient of f at x_k
%       delta: Trust region radius
%   Output:
%       pC: Cauchy point for the model
%   
    % Use pk = g/norm(g) as a first point
    g_norm = norm(g);
    p = - g/g_norm;
    pHp = p'*H*p;
    
    alpha = delta;
    
    % If pBp <= 0, then the minimum is at the edge of the trust region
    % Otherwise, check if the minimum is located inside the trust region
    if pHp > 0
        alpha = min(delta, g_norm/pHp);
    end
    
    pC = 0.99 * alpha * p;
end