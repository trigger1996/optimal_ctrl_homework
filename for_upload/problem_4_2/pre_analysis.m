clear
clc

global A B C Q R E V W
A = [  0,  1;
      -2, -4 ];
B = [ 0 
      0.5 ];
C = eye(2);     % 2 * 2��λ��
Q = [ 4, 0;
      0, 6 ];
R = 0.02;

z = [ 1; 0 ];

E = B * inv(R) * transpose(B);
V = transpose(C) * Q * C;
W = transpose(C) * Q;

Wz = W*z;

% �������ȫ������⣬����ֱ����example 4-1�������ͳ�ʼ״̬�������㣬�õ�����ؽ��
% ϵͳ��PI����4.2��
