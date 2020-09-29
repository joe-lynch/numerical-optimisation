function [xd,md] = dogleg(xn,fx,gx,H,H_,delta)
    
%     This was tested with the brockrosen function, by comparing the results
%     from using the dogleg point vs cauchy point. The results below were
%     with called by
%     tr_dogleg(@f_Rosen,@g_Rosen,eye(2),inv(eye(2)),x0,0.2,1,0.125,1.0e-5),
%     from within the visual function. eta = 10^(-6) in tr_dogleg
%
%     For testing this function seperately, tests such as 
%     [xd,md] = dogleg([1;-1],40,[80;-40], [1,0;0,1], [1,0;0,1], 0.2) can
%     be performed and compared with that of the cauchy value.
%
%     otherwise, the results here similar to the tests in tr_dogleg.m, 
%     x0 = [1,-1]
%     - dogleg : 26 iterations
%     - cauchy : 1246 iterations
%     x0 = [-1.2;1]
%     - dogleg : 47 iterations.
%     - cauchy : 1183 iterations.
%     x0 = [-1.2;1.5]
%     - dogleg : 60 iterations.
%     - cauchy : 1241 iterations.
%     x0 = [-1.3,-0.5]
%     - dogleg : 46 iterations
%     - cauchy : 1272 iterations
%

    curv = gx'*H*gx;
    nrm = norm(gx');
    % Compute cauchy pt and value of model function at that point
    if curv > nrm^3/delta
        alpha = nrm^2/curv;
    else
        alpha = delta/nrm;
    end
    % Set the Cauchy point and evaluate the model at the Cauchy point
    xc = xn  - alpha * gx;
    mc = fx - alpha * nrm^2 + alpha^2/2 * curv;

    % Evaluate newton step
    sN = - H_ * gx;

    if curv <= 0
        % If no curvature info then use cauchy
        xd = xc;
        md = mc;
    else
        if norm(sN) <= delta
        % newton point is inside the trust region, so set dogleg point and
        % model value
        xd = xn + sN;
        md =  - gx'*sN - 0.5*sN'*H*sN;
        else
            % Evaluate unidirectional minimiser step
            sU = -(gx*gx')/curv*gx;
            
            if norm(sU) >= delta
              % the unidirectional minimiser is outside the trust region ( as
              % well as newton point), so use cauchy point
                xd = xc;
                md = mc;
            else
                % unidirectional minimiser is inside, newton point is
                % outside, so find the intersection of the line between
                % them with the trust region radius
                dist = sN - sU;
                a = dist'*dist;
                b = 2*dist'*sU;
                c = sU'*sU - delta^2;
                t  = (-b+sqrt(b^2-4*a*c))/(2*a);
                
                % set the dogleg point and the value of the model
                xd = xn + sU + t*dist;        
                sD = sU+t*dist;
                md =  - gx'*sD - 0.5*sD'*H*sD;
            end
        end
    end
end