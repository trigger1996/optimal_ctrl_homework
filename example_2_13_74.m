clear
clc
clear all

S = dsolve('Dx1 = x2', 'Dx2 = -lambda2', 'Dlambda1 = 0', 'Dlambda2 = -lambda1, x1(0) = 1, x2(0) = 2, x1(2) = 0, lambda2(2) = 0')

t = 0 : 0.01 : 2;
y_x1 = eval(S.x1);
y_x2 = eval(S.x2);
y_u  = -eval(S.lambda2);

plot(t, y_x1, 'r')
hold on
plot(t, y_x2, 'g')
hold on
plot(t, y_u,  'b')