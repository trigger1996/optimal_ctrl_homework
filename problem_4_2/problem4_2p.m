function dp = problem4_2p(t,p)
% Function for P
%
%Define variables to use in external functions
%
%Definition of differential equations
%
% dp=[250*p(2)^2+4*p(2)-2
%     250*p(2)*p(3)-p(1)+3*p(2)+2*p(3)
%     250*p(3)^2-2*p(2)+6*p(3)];
dp = [ 12.5*p(2)^2 + 4*p(2) - 4;
       12.5*p(1)*p(2) - p(1) + 4*p(2) + 2*p(3);
       12.5*p(3)^2 - 2*p(2) + 8*p(3) - 6 ];
