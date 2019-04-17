clear
clc
clear all
x = dsolve('4*x - 2*D2x = 0', 'x(0) = 0, x(2) = 5')
t = 0 : 0.01 : 2;
y = eval(x);        % eval是吧字符串当命令来执行的, 所以必须要有自变量
plot(t,y,'k-')
