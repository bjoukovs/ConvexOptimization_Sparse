close all;clear all; clc

maxit=1e5;
runSubgrad(maxit,0.1);
% hold on
% runSubgrad(maxit,0.5);
% hold on
% runSubgrad(maxit,1);
legend('1/k','10/k','100/k');
title(['Convergence with various stepsizes',' ','MAXIT' ' =',' ',num2str(maxit)])
