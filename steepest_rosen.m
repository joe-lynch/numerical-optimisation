
% Script to apply steepest descent with backtracking line search 
% to the Rosenbrock function

% Set initial guess

x0 = [-1.2;1];
%x0 = [-1.2;1.5];

% Apply basic steepest descent algorithm

theta = 0.1;
tol   = 1.0e-5;

[x,n] = steepest(@f_Rosen,@g_Rosen,x0,theta,tol);

% Set exact solution

xex = [1;1];

% Visualise the solution path 

visual(@f_Rosen,x,x0,xex);

disp('Number of iterations:');
n

for i=1:n+1
    e(i) = norm(x(:,i)-xex);
end

%e(2:n+1)./e(1:n)

disp('Final error:');
e(n+1)

