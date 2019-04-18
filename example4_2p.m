function dp = example4_2p(t,p)
% Function for P
%
%Define variables to use in external functions
%
%Definition of differential equations
%
dp=[25*p(2)^2+4*p(2)-2
    25*p(2)*p(3)-p(1)+3*p(2)+2*p(3)
    25*p(3)^2-2*p(2)+6*p(3)] ;
%%