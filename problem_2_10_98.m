clear
clc
clear all

% dH/du  = u + r2 = 0 ->  u = -r2
% dH/dx1 = -dr1        -> -dr1 = 0
% dH/dx2 = -dr2        -> -dr2 = u
% dH/dr1 =  dx1        ->  dx1 = 0
% dH/dr2 =  dx2        ->  dx2 = u
S = dsolve('-Dr1 = 0', '-Dr2 = -r2', 'Dx2 = -r2', 'Dx1 = x2' , 'x1(0) = 2', 'x1(5) = 0', 'x2(0) = 2')    % 'x2(0) = 2', 'x2(5) = 0'
                                                  % Dx1 = x2 还是Dx1 = 0
                                                  
t = 0 : 0.01 : 5;
y_x1 = eval(S.x1);
y_x2 = eval(S.x2);
y_u  = -eval(S.r2);    % 注意这边有个负号

plot(t, y_x1, 'r')
hold on
plot(t, y_x2, 'g')
hold on
plot(t, y_u,  'b')