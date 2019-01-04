close all;clear all; clc

%Reduced subgradient plots
maxit=200;
%Second argument is stepsize n/it
stepsizes = [0.1 0.5 1 2 4 6];
legends = {};

for i=1:length(stepsizes)

runSubgradReduced(maxit,stepsizes(i)); 
hold on

legends{i} = strcat(num2str(stepsizes(i)),'/k');

end

legend(legends);
title(['Convergence with various stepsizes',' ','MAXIT' ' =',' ',num2str(maxit)])

%%
%Projected subgradient plots
maxit=1e3;
%Second argument is stepsize n/it
runSubgradProj(maxit,1); 
hold on
runSubgradProj(maxit,2);
hold on
runSubgradProj(maxit,3);
legend('1/k','2/k','4/k');
title(['Convergence with various stepsizes',' ','MAXIT' ' =',' ',num2str(maxit)])
