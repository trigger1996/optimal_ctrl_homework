clear
clc
clear all
x = dsolve('4*x - 2*D2x = 0', 'x(0) = 0, x(2) = 5')
t = 0 : 0.01 : 2;
y = eval(x);        % eval�ǰ��ַ�����������ִ�е�, ���Ա���Ҫ���Ա���
plot(t,y,'k-')
