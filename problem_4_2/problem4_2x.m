function dx = problem4_2x(t,x)
% Function for x
%
%Define variables to use in external functions
global tp p tg g
%
%Definition of differential equations
%
%dx=[x(2)
%    -2*x(1)-3*x(2)-250*(interp1(tp,p(:,2),t)*x(1)+ ...
%    interp1(tp,p(:,3),t)*x(2)-interp1(tg,g(:,2),t))];
dx = [ x(2)
       (-12.5 * interp1(tp,p(:,2),t) - 2)*x(1) + (-12.5 * interp1(tp,p(:,3),t) - 4)*x(2) - 12.5*interp1(tg,g(:,2),t) ];
% dx = [ dx1, dx2  ]';
% dx = [a - e*p]*x + e*g
% dx1 =  x2
% dx2 =  (-12.5*p12 - 2) * x1 +(-12.5*p22 - 4) * x2 - 12.5*g2