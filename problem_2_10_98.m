clear
clc
clear all

% 方程初始条件
% H = 0.5u^2 + r1x2 + r2(-2*x1+3*u)
% dH/du  = u + 3r2 = 0          ->  u = -3*r2
% dH/dx1 = -dr1                 -> -dr1 = -2*r2
% dH/dx2 = -dr2                 -> -dr2 = r1
% dH/dr1 =  dx1                 ->  dx1 = x2
% dH/dr2 =  dx2                 ->  dx2 = -2*x1 + 3u

% 根据末端项求解
% S = 0.5 * pi/2 * x1^2
% tf = pi/2
% r1(tf) = dS/dx1               -> pi/2*x1(tf) = r1(tf)
% r2(tf) = dS/dx2               -> 0           = r2(tf)
S = dsolve('-Dr1 = -2*r2', '-Dr2 = r1', 'Dx1 = x2', 'Dx2 = -2*x1 + 3*u', 'x1(0) = 0', 'x2(0)=1', 'r2(pi/2)=0', 'pi/2*x1(pi/2) = r1(pi/2)')

t = 0 : 0.05 : pi/2;
u = -3*S.r2;


% y_x1 = eval(S.x1);
% y_x2 = eval(S.x2);
% y_u  = eval(u);    % 注意这边有个负号
% 
% %plot(t, y_x1, 'r')
% % hold on
% % plot(t, y_x2, 'g')
% % hold on
% % plot(t, y_u,  'b')