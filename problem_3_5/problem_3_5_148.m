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
% https://wenku.baidu.com/view/f40825710b4c2e3f5727635b.html
[ Y, X, t ] = initial(A - B*K, B, C, D, X0);  % initial(A - B*K, BIN, C, D, X0, tfinal);
x1t = [ 1 0 ] * X';
x2t = [ 0 1 ] * X';
ut = -K * X';

figure
plot(t, X(:,1), 'r', 'LineWidth',3)
hold on
plot(t, X(:,2), 'g', 'LineWidth',3)
xlabel('t')
ylabel('x(t)')
l1 = legend('x_1','x_2')
set(l1, 'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',35)
set(gca,'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',25)
 
figure
plot(t, ut, 'b', 'LineWidth',3)
xlabel('t')
ylabel('u(t)')
l2 = legend('Optimal Control')
set(l2, 'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',35)
set(gca,'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',25)
