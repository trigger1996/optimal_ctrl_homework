%%%%%%%%%%%%%
function [x,u]=lqrdnss(As,Bs,Fs,Qs,Rs,x0,kspan)
%
%This m-file calculates and plots the outputs for a
% discrete Linear Quadratic Regulator system
%Based on provided linear state space matrices
% for A and B and Performance Index matrices
% for F, Q and R.
%This function takes these inputs, and using the
% analytical solution to the matrix Riccati equation,
% formulates the optimal states and inputs.
%
%
% SYNTAX: [x,u]=lqrdnss(A,B,F,Q,R,x0,tspan)
%
% INPUTS (All numeric):
% A,B Matrices from xdot=Ax+Bu
% F,Q,R Performance Index Parameters; terminal cost,
% error and control weighting
% x0 State variable initial condition. Must be a
% column vector [xl0;x20;x30 ... ]
% kspan Vector containing sample span [k0 kf]
%
% OUTPUTS:
%  x is the state variable vector
%  u is the input vector
%
% The system plots the Riccati coefficients in
% combinations of 4,
% and the x vector, and u vector in
% combinations of 3.
%
%Check for correct number of inputs
if nargin < 7
    error('Incorrect number of inputs specified')
    return
end

%Convert Variables to normal symbology to
% prevent problems with global statement
A=As;
B=Bs;
F=Fs;
Q=Qs;
R=Rs;
plotflag=0; %set plotflag to a 1 to avoid plotting
% of data on figures
%Define secondary variables for global passing to
% ode-related functions and determine matrice size
[n,m]=size(A);      %Find dimensions of A
[nb,mb]=size(B);    %Find dimensions of B
[nq,mq] =size (Q);  %Find dimensions of Q
[nr,mr]=size(R);    %Find dimensions of R
[nf,mf]=size(F);    %Find dimensions of F
if n~=m             %Verify A is square
    error('A must be square')
else
    [n,n]=size(A);
end
%Data Checks for proper setup
%Check for controllability
if length(A) > rank(ctrb(A,B))
    error('System Not Controllable')
    return
end
if (n ~= nq) || (n ~= mq)
    %Check that A and Q are the same size
    error('A and Q must be the same size');
    return
end
if ~(mf==1 && nf==1)
    if (nq ~= nf) || (mq ~= mf)
        %Check that Q and F are the same size
        error('Q and F must be the same size');
        return
    end
end
if ~(mr==1 && nr==1)
    if (mr ~= nr) || (mb ~= nr)
        error('R must be consistent with B');
        return
    end
end
mq = norm(Q, 1);
% Check if Q is positive semi-definite and symmetric
if any(eig(Q) < -eps*mq) || (norm(Q'-Q,1)/mq > eps)
    disp('Warning: Q is not symmetric and positive semi-definite');
end
mr = norm(R,1);
% Check if R is positive definite and symmetric
if any(eig(R) <= -eps*mr) || (norm(R'-R,1)/mr > eps)
    disp('Warning: R is not symmetric and positive definite');
end
%Define Calculated Matrix
E = B * inv(R) * B';
%Find matrix needed to calculate Analytical Solution
% to Riccati Equation
H = [inv(A), inv(A)*E; Q*inv(A), A' + Q*inv(A)*E];
%Find Eigenvectors
[W, D] =eig(H);
%Find the diagonals from D and pick the negative
% diagonals to create a new matrix M
j=n;
[m1,index1]=sort(real(diag(D)));
for i=1:1:n
    m2(i)=m1(j);
    index2(i)=index1(j);
    index2(i+n)=index1(i+n);
    j=j-1;
end
Md=diag(m2);

%Rearrange W so that it corresponds to the
% sort of the eigenvalues

for i=1:2*n
    w2(:,i)=W(:,index2(i));
end
W=w2;
%Define the Modal Matrix for D and split it into parts
W11=zeros(n);
W12=zeros(n);
W21=zeros(n);
W22=zeros(n);
j=1;
for i=1:2*n:(2*n*n-2*n+1)
    W11(j:j+n-1)=W(i:i+n-1);
    W21(j:j+n-1)=W(i+n:i+2*n-1);
    W12(j:j+n-1)=W(2*n*n+i:2*n*n+i+n-1);
    W22(j:j+n-1)=W(2*n*n+i+n:2*n*n+i+2*n-1);
    j=j+n;
end
% 加一个这个，不然会爆炸
W11 = reshape(W11, n, n);       % transpose(reshape(W11, n, n))
W12 = reshape(W12, n, n);
W21 = reshape(W21, n, n);
W22 = reshape(W22, n, n);

%Find M
M=zeros(n);
j=1;
for i=1:2*n:(2*n*n-2*n+1)
    M(j:j+n-1)=D(i:i+n-1);
    j=j+n;
end

%Zero Vectors
x=zeros(n,1);
%Define Loop Variables (l=lambda)
k0=kspan(1);
kf=kspan(2);
%x and P Conditions
x(:, 1) =x0(:, 1);
Tt=-inv(W22-F*W12)*(W21-F*W11);
P=real((W21+W22*((Md^-(kf-0))*Tt*(Md^-(kf-0)))) * inv(W11+W12*((Md^-(kf-0))*Tt*(Md^-(kf-0)))));
L=inv(R)*B'*(inv(A))'*(P-Q);
u(:,1)=-L*x0(:,1);
k1(1)=0;
for k=(k0+1):1:(kf)
    Tt=-inv(W22-F*W12)*(W21-F*W11);
    P=real((W21+W22*((Md^-(kf-k))*Tt*(Md^-(kf-k)))) * inv(W11+W12*((Md^-(kf-k))*Tt*(Md^-(kf-k)))));
    L=inv(R)*B'*(inv(A))'*(P-Q);
    x(:, k+1)=(A-B*L)*x(:,k);
    u(:, k+1)=-L*x(:,k+1);
    k1(k+1)=k;
end
%Plotting Section, if desired
if plotflag ~= 1
    %Plot Riccati coefficients using flag variables
    % to hold and change colors
    %Variables are plotted one at a time and the plot held
    fig=1; %Figure number
    cflag=1 ;
    j=1;
    Ps=0.;
    for i=1:1:n*n
        %Variable used to change plot color
        %Initialize P Matrix plot variable
        for k=(k0):1:(kf)
            Tt=-inv(W22-F*W12)*(W21-F*W11);
            P=real((W21+W22*(Md^-(kf-k))*Tt*(Md^-(kf-k)))) * inv(W11+W12*((Md^-(kf-k))*Tt*((Md^-(kf-k)))));
            Ps(j)=P(i);
            k2(j)=k;
            j=j+1;
        end
        if cflag==1
            figure(fig);
            plot (k2 , Ps , 'b')
            title('Plot of Riccati Coefficients')
            grid on
            xlabel ('k')
            ylabel('P Matrix')
            hold
            cflag=2;
        elseif cflag==2
            plot(k2 , Ps, 'b')
            cflag=3;
        elseif cflag==3
            plot (k2, Ps, 'b')
            cflag=4;
        elseif cflag==4
            plot (k2 ,Ps, 'b')
            cflag=1;
            fig=fig+1;
        end
    Ps=0.;
    j=1;
    end
if cflag==2||cflag==3||cflag==4
    hold
    fig=fig+1;
end
%Plot Optimized x
x=x' ;
if n>2
    for i=1:3:(3*fix((n-3)/3)+1)
        figure(fig);
        plot(kl,real(x(:,i)),'b',kl,real(x(:,i+l)),'b' ,kl, real(x(:,i+2)),'b');
        grid on
        title('Plot of Optimal States')
        xlabel( 'k')
        ylabel('Optimal States')
        fig=fig+1;
        %
    end
end

if (n-3*fix(n/3))==1
    figure(fig);
    plot(k1,real(x(:,n)),'b')
elseif (n-3*fix(n/3))==2
    figure(fig);
    plot(k1,real(x(:,n-1)),'b',k1,real(x(:,n)),'b')
end
grid on
title('Plot of Optimal States')
xlabel( 'k')
ylabel('Optimal States')
fig=fig+1;
%
%Plot Optimized u
%
u=u' ;
if mb>2
    for i=1:3:(3*fix((mb-3)/3)+1)
    figure(fig);
    plot(k1,real(u(:,i)),'b',k1,real(u(:,i+l)), 'm:',k1,real(u(:,i+2)),'g-.');
    grid on
    title('Plot of Optimal Control')
    xlabel('k')
    ylabel('Optimal Control')
    fig=fig+1;
    end
end
if (mb-3*fix(mb/3))==1
    figure(fig);
    plot(k1,real(u(:,mb)),'b')
elseif (mb-3*fix(mb/3))==2
    figure(fig);
    plot(k1,real(u(: ,mb-1)), 'b' ,k1,real(u(: ,mb)), 'm:')
end
grid on
title('Plot of Optimal Control')
xlabel ('k')
ylabel('Optimal Control')
gtext ('u')
end
%%%%%%%%