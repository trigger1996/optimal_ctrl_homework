clear
clc

global A B C Q R E V W
A = [  0,  1;
      -2, -4 ];
B = [ 0 
      0.5 ];
C = eye(2);     % 2 * 2单位阵
Q = [ 4, 0;
      0, 6 ];
R = 0.02;

z = [ 1; 0 ];

E = B * inv(R) * transpose(B);
V = transpose(C) * Q * C;
W = transpose(C) * Q;

Wz = W*z;

% 这道题完全不好求解，所以直接用example 4-1的条件和初始状态进行运算，得到的相关结果
% 系统和PI都是4.2的
