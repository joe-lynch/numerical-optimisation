
% Script to apply the basic trust region method (using the Cauchy point as
% the approximate solution of the TR subproblem and exact Hessians) to 
% the Rosenbrock function

% Set initial guess

x0 = [1;-1];
%x0 = [-1.2;1.5];

% Apply the basic trust region method with Cauchy point and exact Hessians

tol    = 1.0e-5;
delta0 = 0.5;
delmax = 4;
rho_ac = 0.125;

[x,n] = tr_basic(@f_Rosen,@g_Rosen,@h_Rosen,x0,delta0,delmax,rho_ac,tol);
                       
% Set exact solution

xex = [1;1];

% Visualise the solution path 

visual(@f_Rosen,x,x0,xex);

fprintf('Number of iterations: %i\n',n);

for i=1:n+1
    e(i) = norm(x(:,i)-xex);
end

%e(2:n+1)./e(1:n)

fprintf('Final error: %i\n',e(n+1));
