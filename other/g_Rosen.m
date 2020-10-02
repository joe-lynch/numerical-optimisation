function G = g_Rosen(x)

% The gradient of the two-dimensional Rosenbrock function 

b=10;

G(1) = 2*(x(1)-1) + 4*b*x(1)*(x(1)^2 - x(2));
G(2) = 2*b*(x(2)-x(1)^2);

G = G';

end

