
function [x,n] = tr_basic(f,df,d2f,xn,delta,delmax,rho_ac,tol)

% Iterate up to a maximum of 10000 iterations
    for n=0:10000

        x(:,n+1) = xn;

        % Evaluate the Gradient at the current point x
        gx = df(xn);

        % For the stopping criterion check whether the norm of the gradient is 
        % below the requested tolerance

        if norm(gx) < tol
            break
        end

        % Evaluate objective function and Hessian at the current point x 
        fx = f(xn);
        H  = d2f(xn);
        
        %[xc,mc] = dogleg(xn,fx,gx,H,delta);
        [xc,mc] = cauchy(xn,fx,gx,H,delta);
        
        % Evaluate the model misfit rho_n
        % This is used to see how accurate the new point is

        %  (obj function at xn) - (obj function at xn+1)
        % -----------------------------------------------
        %  (obj function at xn) - (model function at xn+1)

        % The closer it is to 1, the more accurate that the model function was
        % at xn+1

        rho = (fx - f(xc))/(fx - mc);

        % Decide whether to accept or reject the proposed new iterate
        % rho_ac is the predefined minimum accuracy that we require, should be
        % above 0.
        % If true, then set the current point to the new point
        if rho >= rho_ac
            xn = xc;
        end

        %{
          Adjust the trust region radius 
          finds the difference between xn (potentially new point, if it was 
          set), and the current point x(:,n+1)
          Then it finds the norm, to get a distance
          This is the step that is made from the current point to the new point
        %}

        step = norm(xn-x(:,n+1));

        % If the rho that was previously calculated is less than 0.25, then 
        % reduce the size of the tr radius by a quarter, this is because the
        % accuracy of the model at that tr radius was not accurate enough
        if rho < 0.25
            delta = 0.25*delta;
        % Otherwise, if the rho is more than 0.75 (then it was pretty great and
        % accurate) & step size and radius are basically the same, then set the
        % tr radius to double, or the max radius. clearly we need the step size 
        % and radius to be the same, otherwise
        % maybe we found a new point close to the original point, here it
        % wouldn't make much sense to expand the tr.
        elseif (rho > 0.75) && (abs(step-delta) < 1.0e-10)
            delta = min(2*delta,delmax);
        end

        % If the rho was between 0.25 and 0.75, then just leave the tr radius
        % the same size, it was average, neither great nor bad.

    end

end