clear
clc

x1(1) = 1;
x2(1) = 2;

A = [  0,  1
      -1, 1 ];
B = [ 0; 1 ];
Q = [ 1, 0
      0, 0 ];
R = 1;
F = [ 0, 0;
      0, 0 ];

E = B * inv(R) * B';
I = eye(2);

k0 = 0;
kf = 3 + 1;

P_kplus1 = F;
P11(kf) = F(1);
P12(kf) = F(2);
P21(kf) = F(3);
P22(kf) = F(4);
for k = kf - 1 : -1 : 1
    P_k = transpose(A) * P_kplus1 * inv(I + E * P_kplus1) * A + Q;
    P11(k) = P_k(1);
    P12(k) = P_k(2);
    P21(k) = P_k(3);
    P22(k) = P_k(4);
    P_kplus1 = P_k;
end

for k = kf : -1 : 1
    P_k = [ P11(k), P12(k); P21(k), P22(k) ];
    L_k = inv(R) * transpose(B) * inv(transpose(A)) * (P_k - Q);
    l1(k) = L_k(1);
    l2(k) = L_k(2);
end

for k = 1 : kf - 1
    L_k = [ l1(k), l2(k) ];
    x_k = [ x1(k); x2(k) ];
    x_kplus1 = (A - B * L_k) * x_k;
    x1(k + 1) = x_kplus1(1);
    x2(k + 1) = x_kplus1(2);
end

for k = 1 : kf
    L_k = [ l1(k), l2(k) ];
    x_k = [ x1(k); x2(k) ];
    u(k) = -L_k * x_k;
end

%
% plot various values: P(k), x(k), u(k)
% let us first reorder the values of k = 0 to 10
k = 1 : kf;          % 补一句，不然画图会出错
figure(1)
plot(k,P11,'r:o',k,P12,'g:+',k,P22,'b:*', 'linewidth', 3)
xlabel ('k')
ylabel('Riccati Coefficients')
l2 = legend('p_{11}', 'p_{12}', 'p_{22}')
set(l2, 'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',35)
set(gca,'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',25)
% gtext('p_{11}(k)')
% gtext('p_{12}(k)=p_{21}(k)')
% gtext('p_{22}(k)')
%
figure(2)
plot(k,x1,'r:o',k,x2,'g:+', 'linewidth', 3)
xlabel( 'k')
ylabel('Optimal States')
l2 = legend('x_1', 'x_2')
set(l2, 'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',35)
set(gca,'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',25)
% gtext ('x_l (k)')
% gtext ('x_2 (k)')
%
figure(3)
plot (k, u, 'b: *', 'linewidth', 3)
xlabel ('k')
ylabel('Optimal Control')
l3 = legend('u')
set(l3, 'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',35)
set(gca,'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',25)
% gtext('u(k) ')

%%
% t -> infinite
kspan = 20

[Y,L,G]=dlqr(A,B,Q,R);
X(1, :) = [ x1(1), x2(1) ];

for k = 1 : kspan - 1
    X(k + 1, :) = (A - B' * L) * transpose(X(k, :));
end

u_ = []
for k = 1 : kspan
    u_(k,:) = -L * transpose(X(k, :));
end

figure()
plot(1:kspan, X, 'LineWidth', 3)
xlabel ('k')
ylabel('Optimal Control')
l3 = legend('x1', 'x2')
set(l3, 'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',35)
set(gca,'Fontname', 'Times New Roman', 'FontAngle','Italic', 'FontSize',25)

% kspan = 0 : 40;
% x0 = [1; 2];
% [x, u] = lqrdnss(A, B, F, Q, R, x0, kspan);
% 
