clear all, close all;

load('cs.mat');
x_sol = x;

%--------------------------------------
% min norm(x,1) 
% s.t. F_us*x - X_us = 0
%      x >= 0
%--------------------------------------

func = @(x) norm(x,1)

%PROJ = F_us'*inv(F_us*F_us');
PROJ = pinv(F_us);

proj = @(z) z - PROJ*(F_us*z - X_us);

%sub gradient of norm(x,1) = sign(x).
grad = @(z) sign(z);

% Constraint (Ax - b >= 0)
const = @(z) z;
grad_const = @(z) ones(length(z),1);

%%

f = []

MAXIT=10000;
it = 0;

xk = zeros(128,1);

while(it < MAXIT)
    it = it+1;
    
    %Finding a subgradient of norm 1. if the gradient is 0, meaning that
    %x=0, then the subgradient is chosen to be 0.5 to enforce positivity
    
    %The computation of the gradient depends on the inequality constraint
    [cmax, ind] = max(const(xk));
    
    alpha = 0.5/it;
    
    if(cmax <= 0)
        %In this case, the constraint is ok. We use the normal proj subgrad
        %method
        g = grad(xk);
        g = g + (zeros(128,1) == g)/2; %Putting zeros to 0.5
        
    else
        %In this case, the constraint is not ok, xk is not feasible. we set
        %the step direction as the derivative of the constraint
        g = grad_const(xk);
        
    end
    
    %Update
        xk = xk - alpha*proj(g);
    
    
    f = [f func(xk)];
    
end

plot(real(xk))
figure
plot(x_sol)
