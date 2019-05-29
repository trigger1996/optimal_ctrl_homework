clear
clc

% Solution Using Control System Toolbox (STB)
% MATLAB Version 6
%
A=[0.8 1;0,0.6]; %% system matrix A
B=[1;0.5]; %% system matrix B
Q=[1 0;0 1]; %% performance index state weighting matrix Q
R=[1]; %% performance index control weighting matrix R
F=[2,0;0,4]; %% performance index weighting matrix F
%
x1(1)=5; %% initial condition on state xl
x2(1)=3; %% initial condition on state x2
xk=[x1(1);x2(1)];
% note that if kf = 10 then k = [kO,kf] = [0 1 2, ... ,10],
% then we have 11 points and an array xl should have subscript
% x1(N) with N=1 to 11. This is because x(o) is illegal in array

% definition in MATLAB. Let us use N = kf+1
k0=0; % the initial instant k_O
kf=10; % the final instant k_f
N=kf+1; %
[n,n]=size(A); % fixing the order of the system matrix A
I=eye(n); % identity matrix I
E=B*inv(R)*B'; % the matrix E = BR-{-1}B'
%
% solve matrix difference Riccati equation backwards
% starting from kf to kO
% use the form P(k) = A'P(k+1)[I + EP(k+1)]-{-1}A + Q
% first fix the final condition P(k_f) = F
% note that P, Q, R are all symmatric ij = ji
Pkplus1=F;
p11(N)=F(1) ;
p12(N)=F(2);
p21(N)=F(3);
p22(N)=F(4);
%
for k=N-1:-1:1
    Pk = A' *Pkplus1*inv(I+E*Pkplus1)*A+Q;
    p11(k) = Pk(1);
    p12(k) = Pk(2);
    p21(k) = Pk(3);
    p22(k) = Pk(4);
    Pkplus1 = Pk;
end
%
% calcuate the feedback coefficient L
% L = R-{-1}B'A-{-T}[P(k) - Q]
%
for k = N:-1:1
    Pk=[p11(k),p12(k);p21(k),p22(k)] ;
    Lk = inv(R)*B'*inv(A')*(Pk-Q);
    l1(k) = Lk(1);
    l2(k) = Lk(2);
end
%
% solve the optimal states
% x(k+1) = [A-B*L)x(k) given x(O)
%
for k=1:N-1
    Lk = [l1(k),l2(k)];
    xk = [x1(k);x2(k)];
    xkplus1 = (A-B*Lk)*xk;
    x1(k+1) = xkplus1(1) ;
    x2(k+1) = xkplus1(2);
end
%
% solve for optimal control u(k)
%
for k=1:N
    Lk = [l1(k),l2(k)];
    xk = [x1(k);x2(k)];
    u(k) = - Lk*xk;
end
%
% plot various values: P(k), x(k), u(k)
% let us first reorder the values of k = 0 to 10
k = 1 : N;          % 补一句，不然画图会出错
figure(1)
plot(k,p11,'k:o',k,p12,'k:+',k,p22,'k:*')
xlabel ('k')
ylabel('Riccati Coefficients')
gtext('p_{11}(k)')
gtext('p_{12}(k)=p_{21}(k)')
gtext('p_{22}(k)')
%
figure(2)
plot(k,x1,'k:o',k,x2,'k:+')
xlabel( 'k')
ylabel('Optimal States')
gtext ('x_l (k)')
gtext ('x_2 (k)')
%
figure(3)
plot (k, u, 'k: *')
xlabel ('k')
ylabel('Optimal Control')
gtext('u(k) ')
% end of the program
% *********************************************************    