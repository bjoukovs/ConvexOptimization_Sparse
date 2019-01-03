close all;clear all; clc

maxit=1e4;
runSubgrad(maxit,1);
hold on
runSubgrad(maxit,2);
hold on
runSubgrad(maxit,3);
legend('1/k','2/k','4/k');
title(['Convergence with various stepsizes',' ','MAXIT' ' =',' ',num2str(maxit)])
