function [x,n] = tr_dogleg(f,df,H,H_,xn,delta,delmax,rho_ac,tol)
%
%     [x,n] = tr_dogleg(@f_Rosen,@g_Rosen,eye(2),inv(eye(2)),[1;-1],0.2,1,0.125,1.0e-5)
%     x0 = [1;-1]
%     - dogleg + sr1 update    : 26 iterations.
%     - cauchy + sr1 update    : 1246 iterations.
%
%     [x,n] = tr_dogleg(@f_Rosen,@g_Rosen,eye(2),inv(eye(2)),[-0.5;-1],0.2,1,0.125,1.0e-5)
%     x0 = [-0.5;-1]
%     - dogleg + sr1 update    : 41 iterations.
%     - cauchy + sr1 update    : 1161 iterations.
%
%     [x,n] = tr_dogleg(@f_Rosen,@g_Rosen,eye(2),inv(eye(2)),[-1.2;1.5],0.2,1,0.125,1.0e-5)
%     x0 = [-1.2;1.5]
%     - dogleg + sr1 update    : 60 iterations.
%     - cauchy + sr1 update    : 1241 iterations.
%
%     [x,n] = tr_dogleg(@f_Rosen,@g_Rosen,eye(2),inv(eye(2)),[-1.3;-0.5],0.2,1,0.125,1.0e-5)
%     x0 = [-1.3,-0.5]
%     - dogleg + sr1 update    : 46 iterations
%     - cauchy + sr1 update    : 1272 iterations

    for n=0:10000
        x(:,n+1) = xn; 

        % Evaluate the Gradient at the current point x
        gx = df(xn);

        % For the stopping criterion check whether the norm of the gradient is 
        % below the requested tolerance
        if norm(gx) < tol
          break
        end

        % Evaluate objective function
        fx = f(xn);
        
        % Find dogleg point
        [xd,md] = dogleg(xn,fx,gx,H,H_,delta);
        
        % uncomment this, and comment dogleg, to see how it performs with
        % cauchy instead of dogleg...
        %[xd,md] = cauchy(xn,fx,gx,H,delta);

        % Find change in x and gradient values.
        dx = xd - xn; 
        dg = df(xd) - gx; 

        % reduction ratio.
        rho = (fx - f(xd)) / (fx - md);

        if rho >= rho_ac
            xn = xd;
        end

        step = norm(xn-x(:,n+1));
        
        % Update the trust region radius:
        if rho <0.25
            delta = 0.25*delta;
        elseif (rho >0.75) && (abs(step-delta) < 1.0e-10)
            delta = min(2*delta,delmax);
        end
        
        % Small value greater than zero
        eta = 10^(-6);
        
        % Update Hessian and inverse Hessian with SR1
        [H,H_] = sr1(H,H_,dx,dg,eta);
    end
end