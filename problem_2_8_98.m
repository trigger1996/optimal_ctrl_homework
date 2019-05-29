clear
clc
clear all

% 方程初始条件
% H = 0.5*(x1^2 + u^2) + r1*x2 + r2(-2*x1 + 5*u)
% dH/du  = u + 5*r2 = 0         ->  u   = -5*r2
% dH/dx1 = -dr1                 -> -dr1 = x1 - 2*r2
% dH/dx2 = -dr2                 -> -dr2 = r1
% dH/dr1 =  dx1                 ->  dx1 = x2
% dH/dr2 =  dx2                 ->  dx2 = -2*x1 + 5*u
S = dsolve('-Dr1 = x1 - 2*r2', '-Dr2 = r1', 'Dx1 = x2', 'Dx2 = -2*x1 + 5*(-5*r2)', 'x1(0) = 3', 'x1(2) = 0', 'x2(0) = 5', 'x2(2) = 0')