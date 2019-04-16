clear
clc
clear all

A = [ 0., 1.; -2., 1. ];
B = [ 0.; 1. ];
Q = [ 2., 3.; 3., 5. ];
R = [ .25 ];
F = [ 1., 0.5; 0.5, 2 ];
tspan = [ 0  5 ];
x0 = [ 2., -3. ];
[x, u, k] = lqrnss(A, B, F, Q, R, x0, tspan);