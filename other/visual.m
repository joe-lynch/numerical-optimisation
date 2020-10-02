function visual(f,x,x0,xex)
% Postprocessing: Plot contour lines for the objective function f and the 
% solution path of a particular sequence of iterates.

% Specify the window in parameter space that we want to plot.

x1 = [-1.5:0.01:1.5];
x2 = [-1:0.01:1.5];

% Plot the contours of the objective function

Z=zeros(length(x2),length(x1));
for i=1:length(x1)
    for j=1:length(x2)
        Z(j,i) = f([x1(i);x2(j)]);
    end
end
contour(x1,x2,Z,50);
hold on

% Then plot the initial guess, the exact solution and the iterates

plot(x0(1),x0(2),'md');
plot(xex(1),xex(2),'b*');
plot(x(1,:),x(2,:),'r-');

hold off

end

