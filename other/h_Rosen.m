function H = h_Rosen(x)

% The Hessian of the two-dimensional Rosenbrock function 

b = 10;

H(1,1) = 2+12*b*x(1)^2-4*b*x(2);
H(1,2) = -4*b*x(1);
H(2,1) = H(1,2);
H(2,2) = 2*b;

end

