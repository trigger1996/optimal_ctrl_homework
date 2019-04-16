%%%%%%%%%%%%%
%% The following is lqrnss.m
function [x,u,K]=lqrnss(As,Bs,Fs,Qs,Rs,x0,tspan)
%Revision Date 11/14/01
%%
% This m-file calculates and plots the outputs for a
% Linear Quadratic Regulator (LQR) system based on given
% state space matrices A and B and performance index
% matrices F, Q and R. This function takes these inputs,
% and using the analytical solution to the
%% matrix Riccati equation,
% and then computing optimal states and controls.
%
%
%   SYNTAX: [x,u,K]=lqrnss(A,B,F,Q,R,x0,tspan)
%
%   INPUTS (All numeric):
%   A,B            Matrices from xdot=Ax+Bu
%   F,Q,R          Performance Index Parameters;
%   x0             State variable initial condition
%   tspan          Vector containing time span [to tf]
%
%   OUTPUTS:
%
%   x               is the state variable vector
%   u               is the input vector
%   K               is the steady-state matrix inv(R)*B'*P
%
% The system plots Riccati coefficients, x vector,
% and u vector
%
%Define variables to use in external functions
%
global A E F Md tf W11 W12 W21 W22 n,
%
%Check for correct number of inputs
%
if nargin<7
    error('Incorrect number of inputs specified')
    return
end
%
%Convert Variables to normal symbology to prevent
% problems with global statement
%
A=As;
B=Bs;
F=Fs;
Q=Qs;
R=Rs;
plotflag=0;     %set plotflag to 1 to avoid plotting of
% data on figures
%
%Define secondary variables for global passing to
% ode-related functions and determine matrice size
%
[n,m]=size(A);          %Find dimensions of A
[nb,mb]=size(B);        %Find dimensions of B
[nq,mq] =size (Q) ;     %Find dimensions of Q
[nr,mr]=size(R);        %Find dimensions of R
[nf ,mf] =size (F) ;    %Find dimensions of F
if n~=m                 %Verify A is square
    error('A must be square')
else
    [n, n] =size (A) ;
end
%
%Data Checks for proper setup
if length(A) > rank(ctrb(A,B))
    %Check for controllability
    error('System Not Controllable')
    return
end
if (n ~= nq) I (n - mq)
    %Check that A and Q are the same size
    error('A and Q must be the same size');
    return
end
if ~(mf==1&nf==1)
    if (nq ~= nf) | (mq ~= mf)
        %Check that Q and F are the same size
        error('Q and F must be the same size');
        return
    end
end
if ~(mr==1&nr==1)
    if (mr ~= nr) | (mb ~= nr)
        error('R must be consistent with B');
        return
    end
end
mq = norm(Q,1);
% Check if Q is positive semi-definite and symmetric
if any(eig(Q) < -eps*mq) | (norm(Q'-Q,1)/mq > eps)
    disp('Warning: Q is not symmetric and positive semi-definite');
end
mr = norm(R,1);
% Check if R is positive definite and symmetric
if any(eig(R) <= -eps*mr) | (norm(R'-R,1)/mr > eps)
    disp('Warning: R is not symmetric and positive definite');
end
%
%Define Initial Conditions for
%numerical solution of x states
%
t0=tspan (1) ;
tf=tspan(2);
tspan=[tf t0];
%
%Define Calculated Matrices and Vectors
%
E=B*inv(R)*B' ; %E Matrix E=B*(1/R)*B'
%
%Find Hamiltonian matrix needed to use
% analytical solution to
% matrix Riccati differential equation
%
Z=[A,-E;-Q,-A'];
%
%Find Eigenvectors
%
[W, D] =eig (Z) ;
%
%Find the diagonals from D and pick the
% negative diagonals to create
% a new matrix M
%
j=n;
[m1,index1]=sort(real(diag(D)));
    for i=1:1:n
        m2(i)=m1(j);
        index2(i)=index1(j);
        index2(i+n)=index1(i+n);
        j=j-1;
    end
Md=-diag(m2);
%
%Rearrange W so that it corresponds to the sort
% of the eigenvalues
%
for i=1:2*n
    w2(:,i)=W(:,index2(i));
end
W=w2;
%
%Define the Modal Matrix for D and Split it into Parts
%
W11=zeros(n);
W12=zeros(n);
W21=zeros(n);
W22=zeros(n);
j=1 ;
for i=1:2*n:(2*n*n-2*n+1)
    W11(j:j+n-1)=W(i:i+n-1);
    W21(j:j+n-1)=W(i+n:i+2*n-1);
    W12(j:j+n-1)=W(2*n*n+i:2*n*n+i+n-1);
    W22(j:j+n-1)=W(2*n*n+i+n:2*n*n+i+2*n-1);
    j=j+n;
end
%
%Define other initial conditions for
% calculation of P, g, x and u
%
t1=0.;
tx=0.;      %time array for x
tu=0.;      %time array for u      
x=0.;       %state vector
%
%Calculation of optimized x
%
[tx,x]=ode45('lqrnssf',fliplr(tspan),x0, odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6));
%
%Find u vector
%
j=1;
us=0.; %Initialize computational variable
for i=1:1:mb
    for tua=t0:.1:tf
        Tt=-inv(W22-F*W12)*(W21-F*W11);
        P=(W21+W22*expm(-Md*(tf-tua)) * Tt * expm(-Md*(tf-tua))) * inv(W11 + W12 * expm(-Md * (tf-tua)) * Tt * expm(-Md * (tf - tua)));
        K=inv(R)*B'*P;
        xs=interp1(tx,x,tua);
        us1=real(-K*xs');
        us (j) =us1 (i) ;
        tu(j)=tua;
        j=j+1;
    end
    u ( : , i) =us' ;
    us=0;
    j=1 ;
end
%
%Provide final steady-state K
%
P=W21/W11;
K=real(inv(R)*B'*P);
%
%Plotting Section, if desired
%
if plotflag~=1
%
%Plot diagonal Riccati coefficients using a
% flag variable to hold and change colors
%
fig=1 ;
cflag=1;
j=1;
Ps=0. ;
%Figure number
%Variable used to change plot color
%Initialize P matrix plot variable
for i=1:1:n*n
    j = 1
    for t1a=t0: .1:tf
        Tt=-inv(W22-F*W12)*(W21-F*W11);
        P=real((W21+W22*expm(-Md*(tf-t1a))*Tt*expm(-Md * (tf-t1a))) * inv(W11+W12*expm(-Md*(tf-t1a))*Tt *expm(-Md*(tf-t1a))));
        Ps(j)=P(i);
        t1(j)=t1a;
        j=j+1 ;
    end
    if cflag==1
        figure (fig)
        plot(t1,Ps, 'b')
        title('Plot of Riccati Coefficients')
        xlabel ('t')
        ylabel ('P Matrix')
        hold
        cflag=2;
    else
        if cflag==2
            plot (t1, Ps, 'm: ')
            cflag=3;
        elseif cflag==3
            plot(t1,Ps,'g-.')
            cflag=4;
        elseif cflag==4
            plot(t1,Ps,'r--')
            cflag=1 ;
            fig=fig+1;
        end
        Ps=0. ;
        j=1 ;
    end
    if cflag==2 || cflag==3 || cflag==4
        hold
        fig=fig+1;
    end
    %
    %Plot Optimized x
    %
    if n>2
        for i=1:3:(3*fix((n-3)/3)+1)
            figure(fig);
            plot(tx,real(x(:,i)),'b',tx,real(x(:,i+1)),'m:',tx, real(x(:,i+2)),'g-.')
        end
    end
    title('Plot of Optimized x')
    xlabel ('t')
    ylabel('x vectors')
    fig=fig+1;
    if (n-3*fix(n/3))==1
        figure(fig);
        plot(tx,real(x(:,n)),'b')
    else
        if (n-3*fix(n/3))==2
            figure(fig);
            plot (tx, real (x ( : , n -1) ) , 'b' , tx, real (x ( : , n) ) , 'm:')
        end
    title('Plot of Optimized x')
    xlabel ('t')
    ylabel('x vectors')
    fig=fig+1;
    %
    %Plot Optimized u
    %
        if mb > 2
            for i=1:3:(3*fix((mb-3)/3)+1)
                figure(fig);
                plot(tu,real(u(:,i)),'b',tu,real(u(:,i+1)),'m:', tu,real(u(:,i+2)),'g-.')
                title('Plot of Optimized u')
                xlabel ('t')
                ylabel('u vectors')
                fig=fig+1;        
            end
        end
        if (mb-3*fix(mb/3))==1
            figure(fig);
            plot(tu,real(u(:,mb)),'b')
        elseif (mb-3*fix(mb/3))==2
            figure(fig);
            plot(tu,real(u(:,mb-1)),'b',tu,real(u(:,mb)),'m:')            
        end
        title('Plot of Optimized u')
        xlabel (' t')
        ylabel('u vectors')        
    end
end
%
end
%%
%%%%%%%%%%%%%