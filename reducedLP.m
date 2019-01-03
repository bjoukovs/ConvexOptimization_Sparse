clear all, close all;

load('cs.mat');
x_sol = x;

%F_us = real(F_us);

A = [zeros(128,128), F_us];
b = [zeros(128,1); X_us];

F = null(A);

xp = F_us\X_us;

%%

t = 1e6;

ineq_constr = @(z) [-eye(128) -eye(128); -eye(128) eye(128)]*(F*z + xp);

ineq_constr_grad = ([-eye(128) -eye(128); -eye(128) eye(128)]*F)';

func = @(z) ([eye(128), zeros(128,128)]*(F*z + xp)) - 1/t*sum(log(-constr(z)));

grad = @(z) ([eye(128), zeros(128,128)]*F)' - 1/t*sum(  1./(-ineq_constr(z))  .* ineq_constr_grad  ) %% TO CORRECT

%%

maxit = 1e4;
it = 0;
