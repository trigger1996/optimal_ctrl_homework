clear
clc
clear all

% H = (x1^2 + u^2) + r1*x2 + r2(-2*x1 - 3*x2 + 5*u)
% dH/du  = 2*u + 5*r2 = 0       ->  u = -5/2*r2
% dH/dx1 = -dr1                 -> -dr1 = 2*x1 + -2*r2
% dH/dx2 = -dr2                 -> -dr2 = r1 - 3*r2
% dH/dr1 =  dx1                 ->  dx1 = x2
% dH/dr2 =  dx2                 ->  dx2 = -2*x1 - 3*x2 + 5*u
%%
%S = dsolve('-Dr1 = 2*x1 + -2*r2', '-Dr2 = r1 - 3*r2', 'Dx1 = x2', 'Dx2 = -2*x1 - 3*x2 + 5*(-5/2*r2)', 'x1(0)=3', 'x2(0)=2', 'x1(5) = 0', 'x2(5) = 0')


%t = 0 : 0.05 : 10;

%y_x1 = eval(S.x1);
%y_x2 = eval(S.x2);
%y_u  = eval(u);    % 注意这边有个负号

%plot(t, y_x1, 'r')
% hold on
% plot(t, y_x2, 'g')
% hold on
% plot(t, y_u,  'b')

%%
% 2.9同样可以用LQR求解
A = [ 0,  1;
     -2, -3 ];
B = [ 0; 1 ];
Q = [ 1, 0;
      0, 0 ];
R = 1;

[K, P, EV] = lqr(A, B, Q, R);

% X' = AX + Bu = (A - B*K)*X
A - B*K
% dx1 = x2
% dx2 = -2.2361 * x1 - 3.0777
BIN = [0, 0; 0, 0];
C   = [1,0;0,1];
D   = [1,0;0,1];

X0 = [3, 2];
[ Y, X, t ] = initial(A - B*K, BIN, C, D, X0);
ut = -K * X';

figure
plot(t, X(:,1), 'r')
hold on
plot(t, X(:,2), 'g')

figure
plot(t, ut, 'b')
xlabel('t')
ylabel('u(t)')
