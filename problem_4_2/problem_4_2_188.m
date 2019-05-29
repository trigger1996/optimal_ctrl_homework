clear all

% 计算A, B, Q, R
pre_analysis
%
%Define variables to use in external functions
global tp p tg g
%
%Define Initial Conditions for numerical solution
% of g and x states
%
tf=20;
tspan=[tf 0];
tp=0.;
tg=0.;
tx=0.;
pf=[2. ,0. ,0.];
gf=[2. ,0.] ;
xO=[-0.5,0.] ;
%
%Calculation of P
%
[tp,p]=ode45('problem4_2p',tspan,pf,odeset('refine',2, ...      % tspan是时间，pf是初始值，关于odeset: http://www.ece.northwestern.edu/CFS/local-apps/matlabhelp/techdoc/ref/odeset.html
'RelTol',1e-4,'AbsTol',1e-6));                                  % tp是返回列向量的时间点
%
%Calculation of g
%
[tg,g]=ode45('problem4_2g',tp,gf,odeset('refine',2, ...
'RelTol',1e-4,'AbsTol',1e-6));
%
%Calculation of optimized x
%
[tx,x]=ode45('problem4_2x',flipud(tg),xO, ...
odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6));
%
%Plot Riccati coefficients
%
fig=1; %Figure number
figure (fig)
plot(tp,real(p(:,1)),'r',tp,real(p(:,2)),'g',tp, ...
real (p ( : ,3) ) , 'b')
grid on
xlabel (' t')
ylabel('Riccati Coefficients')
hold
%
fig=fig+1;
%
%Plot g values
%
figure(fig);
plot(tg,real(g(:,1)),'r',tg,real(g(:,2)),'g')
grid on
xlabel (' t')
ylabel('g vector')
%%
%
fig=fig+1;
%
%Plot Optimal States x
%
figure(fig);
plot(tx,real(x(:,1)),'r',tx,real(x(:,2)),'g')
grid on
xlabel (' t')
ylabel('Optimal States')
%
fig=fig+1;
%
%Plot Optimal Control u
%
[n,m] =size (p) ;
p12=p(: ,2) ;
p22=p ( : , 3) ;
x1=x(: ,1) ;
x2=x(: ,2);
g2=flipud(g(:,2));
for i=1:1:n
u(i) = -250*(p12(i)*x1(i) + p22(i)*x2(i) - g2(i));
end
figure(fig);
plot(tp,real(u),'b')
grid on
xlabel ('t')
ylabel('Optimal Control')