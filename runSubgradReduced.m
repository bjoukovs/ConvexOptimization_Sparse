function [numit,x] = runSubgradReduced(MAXIT,step)
load('cs.mat');

%declarations
x_sol = x;
Fopt=3;
idx=0;
iter=[];
it = 0;
xk = zeros(72,1);

%Decomposition of the problem in real and imaginary parts
X_us2 = [real(X_us); imag(X_us); zeros(128,1)];
F_us2 = [real(F_us) -imag(F_us); imag(F_us) real(F_us); zeros(128,128) eye(128)];

%Particular solution
xp = pinv(F_us2)*X_us2;

%F condition satisfied
F = null(F_us2);

%--------------------------------------
% min norm(x,1) 
% s.t. F_us*x - X_us = 0
%      imag(x) = 0
%      real(x) >= 0
%
%--------------------------------------
rmask = [eye(128) zeros(128,128); zeros(128,128) zeros(128,128)];
imask = [zeros(128,128) zeros(128,128); zeros(128,128) eye(128)];

%Get real and imaginary parts
get_x_real = @(z) rmask*(F*z + xp);
get_x_im = @(z) imask*(F*z + xp);

%Objective function
func = @(z) norm(get_x_real(z)+ 1i*get_x_im(z), 1);

%sub gradient of norm(x,1) = sign(x).
grad = @(z) F'*sign(F*z + xp);

% Constraint (Ax - b <= 0)
const = @(z) -get_x_real(z);
grad_const = @(z) -(rmask*F)';
    
F = [];

BEST_XK = xk;

while(it < MAXIT)
    it = it+1;
    
    %Finding a subgradient of norm 1. if the gradient is 0, meaning that
    %x=0, then the subgradient is chosen to be 0.5 to enforce positivity
    
    %The computation of the gradient depends on the inequality constraint
   
    %[cmax, ind] = max(const(xk));
    
    %Finding most violated contraint, should be positive if violated
    [cmax, ind] = max(const(xk));
    
    j = 0;
    
    while(j < 1000 && cmax >= 1e-3)
        gc = grad_const(xk);
        g = gc(:,ind);
       
        stepsize = abs(cmax)/norm(g)^2;
        
        xk = xk-stepsize*g;
        
        [cmax, ind] = max(const(xk));
        j = j+1;
    end
    
    %In this case, the constraint is ok. We use the normal subgrad for the
    %update  
    g = grad(xk);
    
    stepsize = step/it;

    %Update
    xk = xk - stepsize*g; 
    
    F = [F func(xk)];
    
    %Check for best point
    if func(xk) < func(BEST_XK)
        BEST_XK = xk;
    end
    
    objective(it)=func(BEST_XK)-Fopt;
    
    x2 = get_x_real(BEST_XK);
    x2=x2(1:128);
    
    %Calculate RMSE based # of iterations given tol
    RMSE(it) = sqrt(mean((x2 - x_sol).^2));
    
    if RMSE(it) < 1e-3%5.5333e-04
        idx=idx+1;
        iter(idx)=it;
    end
    
end
    
numit=iter(1);
x = get_x_real(BEST_XK);

fopt=0;
semilogy((objective(1:end)-fopt),'LineWidth',1.5)
xlabel('Number of Iterations');ylabel('f(x_k)-f^*');
grid on
 




