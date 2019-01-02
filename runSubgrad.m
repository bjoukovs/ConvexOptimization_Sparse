function [pbjective] = runsubgrd(MAXIT,step)
%clear all, close all;clc

load('cs.mat');
fopt=3;
x_sol = x;

%Decomposition of the problem in real and imaginary parts
X_us2 = [real(X_us); imag(X_us)];
F_us2 = [real(F_us) -imag(F_us); imag(F_us) real(F_us)];

%--------------------------------------
% min norm(x,1) 
% s.t. F_us*x - X_us = 0
%      real(x) >= 0
%
%--------------------------------------

func = @(x) norm(x(1:128) + 1i*x(129:256), 1);

%PROJ = F_us'*inv(F_us*F_us');
PROJ = pinv(F_us2);
proj = @(z) z - PROJ*(F_us2*z - X_us2);

%sub gradient of norm(x,1) = sign(x).
grad = @(z) [sign(z(1:128)); zeros(128,1)];

% Constraint (Ax - b <= 0)
const = @(z) -[z(1:128); zeros(128,1)];
grad_const = @(z) -[ones(128,1); zeros(128,1)];

f = [];
%MAXIT=1e5;
%MAXIT=1000;
it = 0;
xk = zeros(256,1);
%incr=1;
test=0;

    
while(it < MAXIT)
    it = it+1;
    
    %Finding a subgradient of norm 1. if the gradient is 0, meaning that
    %x=0, then the subgradient is chosen to be 0.5 to enforce positivity
    
    %The computation of the gradient depends on the inequality constraint
   
    %[cmax, ind] = max(const(xk));
    
    %Finding most violated contraint, should be positive if violated
    [cmax, ind] = max(const(xk));

    
    if(cmax <= 0)
        %In this case, the constraint is ok. We use the normal proj subgrad
        %method
        g = grad(xk);
        g = g - ([zeros(128,1); 10*ones(128,1)] == g)/2; %Putting zeros to 0.5
        
        switch step
            case 1
                stepsize = 1/it;          
            case 2
                stepsize = 2/(it);                
            case 3
                stepsize = 4/(it); 
                 
        end
        
        %Update
        xk = xk - stepsize*proj(g);  
        
      
        
         
    else
        %In this case, the constraint is not ok, xk is not feasible. we set
        %the step direction as the derivative of the constraint
        g = zeros(256,1);
        gc = grad_const(xk);
        g(ind) = gc(ind);
        
        stepsize = abs(cmax)*2;
        
        xk = xk-stepsize*g;
        
    end
    
    
    objective(it)=func(xk);
    
    
end

fopt=0;
 semilogy((objective(2:end)-fopt),'LineWidth',1.5)
 xlabel('Number of Iterations');ylabel('f(x_k)-f^*');
 grid on
 %hold on





% subplot(2,1,1)
% plot(xk(1:128))
% subplot(2,1,2)
% plot(x_sol)


end

