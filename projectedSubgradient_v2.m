clear all, close all;

load('cs.mat');
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
    
MAXIT=50000;
%MAXIT=100;
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
    

%     if abs(posm) > abs(negm)
%         cmax=posm;
%     else
%         cmax=negm;   
%     end
    

    
    %cmax = max(const(xk));
    %test2(:,it)=(max(abs(real(const(xk)))));

    %Better with the stepsize rule ak=0.1/k
    alpha = 0.5/it;
    
    if(cmax <= 0)
        %In this case, the constraint is ok. We use the normal proj subgrad
        %method
        g = grad(xk);
        g = g - ([zeros(128,1); 10*ones(128,1)] == g)/2; %Putting zeros to 0.5
        
        %Update
        xk = xk - alpha*proj(g);  
        
         
    else
        %In this case, the constraint is not ok, xk is not feasible. we set
        %the step direction as the derivative of the constraint
        g = zeros(256,1);
        gc = grad_const(xk);
        g(ind) = gc(ind);
        
        alpha = abs(cmax)*2;
        
        xk = xk-alpha*g;
        
    end
    
    
end

subplot(2,1,1)
plot(xk(1:128))
subplot(2,1,2)
plot(x_sol)