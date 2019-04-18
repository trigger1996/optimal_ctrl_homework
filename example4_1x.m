function dx = example4_1x(t,x)
% Function for x
%
%Define variables to use in external functions
global tp p tg g
%
%Definition of differential equations
%
dx=[x(2)
    -2*x(1)-3*x(2)-250*(interp1(tp,p(:,2),t)*x(1)+ ...
    interp1(tp,p(:,3),t)*x(2)-interp1(tg,g(:,2),t))];
%