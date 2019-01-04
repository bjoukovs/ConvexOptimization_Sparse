close all;clear all; clc
load('cs.mat');
%Reduced subgradient plots
maxit=1000;
%Second argument is stepsize n/it
%stepsizes = [0.1 0.2 0.5 1];
stepsizes = (0.1:0.1:1);

legends = {};

for i=1:length(stepsizes)

[numit(i),xsolved]=runSubgradReduced(maxit,stepsizes(i)); 
hold on

legends{i} = strcat(num2str(stepsizes(i)),'/k');

legend(legends);
title(['Convergence with various stepsizes',' ','MAXIT' ' =',' ',num2str(maxit)])


end

% figure
% subplot(2,1,1)
% plot(xsolved(1:128));
% legend('Reduced subgradiant')
% subplot(2,1,2);
% plot(x);
% legend('True solution');

%sgtitle('Comparison of solutions')


%%
plot(objective);

%%
plot(stepsizes,numit,'b','LineWidth',2);
ylabel('# Iterations');xlabel('Stepsize \alpha')
grid on
%%
maxit=100;
[numit,objective]=runSubgradReduced(maxit,0.1); 

%%
%Projected subgradient plots
maxit=1000;
%Second argument is stepsize n/it
stepsizes = [0.1 0.2];
%stepsizes = (0.1:0.1:1);

legends = {};

for i=1:length(stepsizes)
[numit(i)]=runSubgradReduced(maxit,stepsizes(i)); 
hold on

legends{i} = strcat(num2str(stepsizes(i)),'/k');

end
