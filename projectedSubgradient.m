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
grad = @(z) sign(real(z));

% Constraint (Ax - b >= 0)
const = @(z) real(z);
grad_const = @(z) ones(length(z),1);

f = [];
    
MAXIT=10000;
%MAXIT=100;
it = 0;
xk = zeros(128,1);
%incr=1;
test=0;
while(it < MAXIT)
    it = it+1;
    
    %Finding a subgradient of norm 1. if the gradient is 0, meaning that
    %x=0, then the subgradient is chosen to be 0.5 to enforce positivity
    
    %The computation of the gradient depends on the inequality constraint
   
    %[cmax, ind] = max(const(xk));
    
    %Finding most violated contraint, should be negative if violated
    cmin = min(const(xk));
    

%     if abs(posm) > abs(negm)
%         cmax=posm;
%     else
%         cmax=negm;   
%     end
    

    
    %cmax = max(const(xk));
    %test2(:,it)=(max(abs(real(const(xk)))));

    %Better with the stepsize rule ak=0.1/k
    alpha = 1/it;
    
    if(cmin >= 0)
        %In this case, the constraint is ok. We use the normal proj subgrad
        %method
        g = grad(xk);
        g = g - (zeros(128,1) == g)/2; %Putting zeros to 0.5
        %test(it,:)=(zeros(128,1) == g)/2;
        
         
    else
        %In this case, the constraint is not ok, xk is not feasible. we set
        %the step direction as the derivative of the constraint
        g = grad_const(xk);
       
        %testblu(it)=1+1;
        %g = grad_const(randi([-1 1]));
        
    end
    
    %Update
    xk = xk - alpha*proj(g);  
end

subplot(2,1,1)
plot((real(xk)))
subplot(2,1,2)
plot(x_sol)