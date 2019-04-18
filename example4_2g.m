function dg = example4_2g(t,g)
% Function for g
%
%Define variables to use in external functions
%
global tp p
%
%Definition of differential equations
%
dg=[(25*interp1(tp,p(:,2),t)+2)*g(2)-4*t
    -g(1)+(25*interp1(tp,p(:,3),t)+3)*g(2)] ;
%%