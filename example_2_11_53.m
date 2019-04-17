clear
clc
clear all
x = dsolve('D2x - 2*x = 0', 'x(0) = 1, x(1) = 0')
t = 0 : 0.01 : 1;
y = eval(x);        % eval�ǰ��ַ�����������ִ�е�, ���Ա���Ҫ���Ա���
plot(t,y,'k-')

%dsolve��ô��ͼ��ʵ��
% y = dsolve('x*D2y - 3*Dy ==x^2','y(1)=0','y(5) == 0','x');
% x=-2*pi:pi/10:2*pi;
% y1=eval(y);
% plot(x,y1,'k-')