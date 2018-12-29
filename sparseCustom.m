clear all, close all;

load('cs.mat');
x_sol = x;


% Variable vector [x v]
x = [zeros(128, 1); zeros(128, 1)];

%% 1 st step: rewrite the problem by eliminating the constraints

% Equality constraint [F_us 0][x v]' = X_us
A = [F_us zeros(128,128)];
b = X_us;

F = null(A);
xp = A\b;

% New variable z
sz = size(F);
z = zeros(sz(2), 1);
z = randn(sz(2), 1);


% Inequality constraints: use log-barrier trick

% defining the unconstrained problem
%---------------------
% min [1 0](F*z + xp) - 1/t*sum(log(-[I I; I -I]*(F*z + xp)))
%---------------------

t = 10e3;
A_constr = [eye(128), eye(128); eye(128), -eye(128)];

cost = @(z) [ones(1,128) zeros(1,128)]*(F*z + xp) - 1/t*sum( log(-A_constr*(F*z + xp)) );

% Gradient dF/dz = ([1 0]F)' - 1/t*

P = @(z) -A_constr*F*z - A_constr*xp;

grad = @(z) ([ones(1,128) zeros(1,128)]*F)' - 1/t* sum( repmat(1./P(z), 1, sz(2)).*(-A_constr*F) , 1)';

 %%

[x, it, prec] = backtracking_gradient_descent(z, 0.5, 0.5, cost, grad);

solution = F*x + xp;
plot(solution);
