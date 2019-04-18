clear all
%
%Define variables to use in external functions
global tp p tg g
%
%Define Initial Conditions for numerical solution of
% g and x states
%
tf=20;
tspan=[tf 0];
tp=0.;
tg=0.;
tx=0.;
pf=[0. ,0. ,0.];
gf = [0. , 0.] ;
xO= [-1. , 0.] ;
%
%Calculation of P
[tp,p]=ode45('example4_2p',tspan,pf, ...
odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6));
%
%Calculation of g
%
[tg,g]=ode45('example4_2g',tp,gf, ...
odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6));
%
%Calculation of optimized x
%
[tx,x]=ode45('example4_2x',flipud(tg),xO, ...
odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6));
%
fig=1; %Figure number
figure (fig)
plot(tp,real(p(:,1)),'k',tp,real(p(:,2)),'k',tp, ...
    real(p(:,3)),'k')
grid on
title('Plot of P')
xlabel ('t')
ylabel('Riccati Coefficients')
hold
%
fig=fig+1;
%
%Plotg values
%
figure(fig);
plot(tg,real(g(:,1)),'k',tg,real(g(:,2)),'k')
grid on
title('Plot of g Vector')
xlabel ('t')
ylabel('g vector')
%
fig=fig+1;
%
%Plot Optimized x
%
figure(fig);
plot(tx,real(x(:,1)),'k',tx,real(x(:,2)),'k')
grid on
title('Plot of Optimal States')
xlabel ('t')
ylabel('Optimal States')
%
fig=fig+1 ;
%
%Calculate and Plot Optimized u
%
[n,m]=size(p);
p12=flipud(p(:,2));
p22=flipud(p(:,3));
x1=x(: ,1) ;
x2=x(: ,2);
g2=flipud(g(:,2));