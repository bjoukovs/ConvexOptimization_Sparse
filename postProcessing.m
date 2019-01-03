close all;clear all; clc

maxit=1e2;
x = runSubgrad(maxit,1);
hold on
x2 = runSubgrad(maxit,2);
hold on
x3 = runSubgrad(maxit,3);
legend('1/k','2/k','4/k');
title(['Convergence with various stepsizes',' ','MAXIT' ' =',' ',num2str(maxit)])

figure
plot(x)
hold on
plot(x2)
plot(x3)