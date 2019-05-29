clear
clc
clear all

A = [  0, 1;
      -2, 4 ];
B = [ 0; 5 ];
Q = [ 5, 0;
      0, 2 ];
R = 4;

[K, P, EV] = lqr(A, B, Q, R)