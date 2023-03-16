clear
x0 = [0.15;0]; % initial point for x
lambda0 = 0; % initial point for lambda
numIter = 1000; % number of iterations
stepSize = 0.01; % stepsize
x = zeros(numIter+1,2);
h=0.01;
k=x-h;
lambda = zeros(numIter+1,1);
x(1,:) = x0; % initial point for x
lambda(1) = lambda0; % initial point for lambda
for i = 2:numIter+1
    f = (4-2.1*x(i-1,1)^2+x(i-1,1)^4/3)*x(i-1,1)^2+x(i-1,1)*x(i-1,2)+(-4+4*x(i-1,2)^2)*x(i-1,2)^2;
    %f=f(x,y)
    fx2 = (4-2.1*k(i-1,1)^2+k(i-1,1)^4/3)*k(i-1,1)^2+k(i-1,1)*k(i-1,2)+(-4+4*k(i-1,2)^2)*k(i-1,2)^2;
    %fx2=f(x-h)
    fx=(f-fx2)/h;
    %fx=(f(x)-f(x-h))/h, finit difference of f(x)
    fy2=(4-2.1*x(i-1,1)^2+x(i-1,1)^4/3)*x(i-1,1)^2+x(i-1,1)*k(i-1,2)+(-4+4*k(i-1,2)^2)*k(i-1,2)^2;
    %fy2=f(y-h)
    fy = (f-fy2)/h;
    %finit difference of f(y)
    hx = 1;
    hy = -1;
    gradf = [fx;fy]; % gradient of f
    gradh = [hx;hy]; % gradient of h
    x(i,:) = x(i-1,:) - stepSize*(gradf' + lambda(i-1)*gradh'); % update x
    lambda(i) = lambda(i-1) + stepSize*(x(i-1,1)-x(i-1,2)-1); % update lambda
end
nonlinearFunctionConstrained; % draw the function
figure(1)
term1 = (4-2.1*x(end,1).^2+x(end,1).^4/3).*x(end,1).^2;
term2 = x(end,1).*x(end,2);
term3 = (-4+4*x(end,2).^2).*x(end,2).^2;
optf = term1 + term2 + term3; % function value at the optimal point
hold on,
plot3(x(end,1),x(end,2),optf+0.2,'ro','LineWidth',8)
