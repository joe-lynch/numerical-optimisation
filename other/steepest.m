
function [x,n] = steepest(f,df,xn,theta,tol)

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
    
    % Evaluate the objective function at the current point x
    
    fx = f(xn);

    % Call backtracking line search with steepest descent direction
    
    alpha = linesearch(xn,fx,gx,-gx,theta,f);
    
    % Calculate the next iterate in steepest descent direction
    
    xn  = xn - alpha*gx;
              
end

end