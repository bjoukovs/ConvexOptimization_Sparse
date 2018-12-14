clear all;
close all;

load('cs.mat');

cvx_begin
    variable xs(128)
    minimize( norm(xs, 1) )
    subject to
        F_us*xs == X_us
cvx_end



% k = 5;
% cvx_begin
%     variable xs1(128)
%     minimize( norm(F_us*xs1 - X_us, 2) )
%     subject to
%        norm(xs1,1) <= k
%        xs1 >= zeros(128,1)
% cvx_end

% gamma = 1e-3;
% 
% cvx_begin
%     variable xs2(128)
%     minimize( norm(F_us*xs2 - X_us, 2) + gamma*norm(xs2, 1))
%     
%     subject to
%         xs2 >= zeros(128,1)
% cvx_end


cvx_begin
    variable v(128)
    variable xs3(128)
    
    minimize( ones(1,128)*v )
    
    subject to
        F_us*xs3 - X_us == 0
        %xs3 >= - v
        %xs3 <= v
        [eye(128) eye(128); -eye(128) eye(128)]*[xs3;v] >= 0
cvx_end


figure
plot(x)
%figure
%plot(xs)

figure
plot(xs3)