function dp = example4_1p(t,p)
% Function for P
%
%Define variables to use in external functions
%
%Definition of differential equations
%
dp=[250*p(2)^2+4*p(2)-2
%
250*p(2)*p(3)-p(1)+3*p(2)+2*p(3)
250*p(3)^2-2*p(2)+6*p(3)];