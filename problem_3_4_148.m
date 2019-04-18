clear
clc
clear all

x10 = 3;
x20 = 5;
X0  = [ x10; x20 ];

A = [  0, 1;
      -2, -3 ];
B = [ 0; 1 ];
Q = [ 1, 0;
      0, 1 ];
R = 1;

[K, P, EV] = lqr(A, B, Q, R)

BIN = [ 0; 0 ];
C   = [ 1  1 ];
D   = [ 1 ];
tfinal = 10;
t = 0 : 0.05 : tfinal;
[ Y, X, t ] = initial(A - B*K, BIN, C, D, X0, tfinal);
x1t = [ 1 0 ] * X';
x2t = [ 0 1 ] * X';
ut = -K * X';
figure;
plot(t, x1t, 'r', t, x2t, 'g');
figure;
plot(t, ut, 'b')
xlabel('t')
ylabel('u(t)')
