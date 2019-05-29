%% the MATLAB. Version 6
%% The following file example.m requires
%% two other files lqrnss.m and lqrnssf.m
%% which are given in Appendix C
clear all
clc
A= [ .8 , 1; 0, .6] ;
B= [1 ; .5] ;
F= [2 , 0 ; 0 ,4] ;
Q=[1,0;0,1] ;
R=1;
kspan=[0 10];
x0( : ,1) = [5. ; 3.] ;
[x,u]=lqrdnss(A,B,F,Q,R,x0,kspan);

% 如果是离散的
%C   = [ 1  1 ];
%D   = [ 1 ];
%[X,L,G]=dare(A,B,Q,R);
%[ Y, X, t ] = initial(A - B*transpose(real(L)), B, C, D, x0);  % initial(A - B*K, BIN, C, D, X0, tfinal);
%figure;
%plot(t, X)