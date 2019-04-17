clear
clc
clear all

X0  = 3;        % 0
A   = 1;
B   = 1;
Q   = 2;
R   = 0.25;

[K, P, EV] = lqr(A, B, Q, R)

BIN = 0;
C   = 1;
D   = 1;
tfinal = 10;
t = 0 : 0.05 : tfinal;
[ Y, X, t ] = initial(A - B*K, BIN, C, D, X0, tfinal);
ut = -K * X;
figure
plot(t, X, 'r')
figure
plot(t, ut, 'b')
xlabel('t')
ylabel('u(t)')