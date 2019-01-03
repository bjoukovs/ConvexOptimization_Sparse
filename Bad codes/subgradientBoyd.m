%%
close all;clear all; clc
load('cs.mat');
% Solve linear program
%   minimize    c'*x
%   subject to  Ax <= b
% using subgradient method
%
% generates plot of best objective value versus iterations
%
% EE364b Convex Optimization II, S. Boyd
% Written by Almir Mutapcic, 01/19/07
%

A = F_us;
b = X_us;     % positive coefficients, so zero is a feasible point
c = -A'*ones(128,1);
%c=ones(128,1);
m=length(A);n=m;

%Optimal value of LP from CVX solution
f_min = 3;  

%********************************************************************
% subgradient method
%********************************************************************
f = [+Inf]; fbest = [+Inf]; fconstr = [];

k = 1;
MAX_ITERS = 20000;
EPS = 1e-3; % we'll overstep this distance in feasibility steps

% initial point
xp = zeros(n,1);

while k < MAX_ITERS 
  % feasibility check
  [fval,ind] = max(A*xp - b);

  % subgradient and step size computation
  if( fval > 0 ) % feasibility step
    fbest(end+1) = fbest(end);
    g = A(ind,:)';
    alpha = (fval+EPS)/norm(g)^2;

  else % optimality step
    f(end+1) = c'*xp;
    fbest(end+1) = min( c'*xp, fbest(end) );
    g = c;
    alpha = 1/k;
  end

  % constraint violation values
  fconstr(end+1) = fval;

  % subgradient update
  xp = xp - alpha*g; k = k + 1;

  % print out fbest every 100 iterations
  if( rem(k,100) == 0 ), fprintf(1,'iter: %d\n',k), end
end

subplot(2,1,1)
plot(x);
subplot(2,1,2)
plot(real(xp))