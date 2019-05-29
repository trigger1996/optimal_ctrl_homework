clear
clc

% H = x^2 + u^2 + r1*(-x + u)
% dH/du = 0                --> 2*u = -r
% dH/dx = -dr              --> -dr = 2*x - r
% dH/dr = dx               -->  dx = -x + u

% ¡Ótf = 5
S = dsolve('-Dr = 2*x - r', 'Dx = -x + (-r/2)', 'x(0) = 5')
SS = solve('x^2 + (r^2)/4 + r*(-x + -r/2) = 0', strcat('x=', string(S.x)), strcat('r=', string(S.r)), 'C3', 't')

% u = -S.r/2;
% t = 0 : 0.05 : 10;
% 
% y_x = eval(S.x);
% y_r = eval(S.r);
% y_u = -y_r / 2;
%  
% plot(t, y_x, 'r')
% hold on
% plot(t, y_r, 'g')
% hold on
% plot(t, y_u, 'b')
% 
