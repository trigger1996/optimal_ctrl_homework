clear
clc
clear all

x10 = 2;
x20 = 3;
X0  = [ x10; x20 ];

A   = [ 0, 1; -2, 1];
B   = [ 0; 1 ];
Q   = [ 2, 3; 3, 5 ];
R   = [ 0.25 ];

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
plot(t, x1t, 'r', t, x2t, 'g');
gtext('x_1(t)')
gtext('x_2(t)')
plot(t, ut, 'b')
xlabel('t')
ylabel('u(t)')
