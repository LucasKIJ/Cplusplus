clear all
close all
data = load('xy.dat');
x = linspace(0,1);
plot(x,exp(-x),'linewidth',2)
hold on
plot(data(:,1),data(:,2),'linewidth',2)
legend('true','Euler')