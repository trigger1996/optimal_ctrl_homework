function dg = problem4_2g(t,g)
% Function for g
%
%Define variables to use in external functions
%
global tp p
%
%Definition of differential equations
%
dg=[(12.5*interp1(tp,p(:,2),t) + 2)*g(2) - 4
    -g(1)+(12.5*interp1(tp,p(:,3),t) + 4) * g(2)];
% dg1 = (12.5 * p12 + 2) * g2 - 4
% dg2 = -g1 + (12.5 * p22 + 4) * g2