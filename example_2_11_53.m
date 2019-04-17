clear
clc
clear all
x = dsolve('D2x - 2*x = 0', 'x(0) = 1, x(1) = 0')
t = 0 : 0.01 : 1;
y = eval(x);        % eval是吧字符串当命令来执行的, 所以必须要有自变量
plot(t,y,'k-')

%dsolve怎么画图的实例
% y = dsolve('x*D2y - 3*Dy ==x^2','y(1)=0','y(5) == 0','x');
% x=-2*pi:pi/10:2*pi;
% y1=eval(y);
% plot(x,y1,'k-')