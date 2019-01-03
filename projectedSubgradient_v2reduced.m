clear all, close all;

load('cs.mat');
x_sol = x;

%Decomposition of the problem in real and imaginary parts
X_us2 = [real(X_us); imag(X_us); zeros(128,1)];
F_us2 = [real(F_us) -imag(F_us); imag(F_us) real(F_us); zeros(128,128) eye(128)];

%Particular solution
xp = pinv(F_us2)*X_us2;

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
get_x_real = @(z) rmask*(F*z + xp);
get_x_im = @(z) imask*(F*z + xp);


func = @(z) norm(get_x_real(z)+ 1i*get_x_im(z), 1);

%PROJ = F_us'*inv(F_us*F_us');
%PROJ = pinv(F_us2);

%sub gradient of norm(x,1) = sign(x).
%grad = @(z) [sign(z(1:128)); zeros(128,1)];
grad = @(z) F'*sign(F*z + xp);

% Constraint (Ax - b <= 0)
const = @(z) -get_x_real(z);
grad_const = @(z) -(rmask*F)';

f = [];
    
MAXIT=1e4;
%MAXIT=100;
it = 0;
xk = zeros(72,1);
%xk = pinv(F_us2)*X_us2;
%incr=1;
test=0;


alpha = 0.5;
beta = 0.1;

J = [];
K = [];

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
        %g = g - ([zeros(128,1); 10*ones(128,1)] == g)/10; %Putting zeros to -0.1 to enforce positivity
        
        stepsize = 1/it;
        
         t = 1;
        j=0;
        
        while(j<10 && func(xk - alpha*t*g) < func(xk) - alpha*t*g'*g)
            
            t = beta*t;
            j = j+1;
            
        end

        
        
        stepsize = alpha*t;
        J = [J stepsize];
        K = [K 1];

        %Update
        xk = xk - stepsize*g; 

         
    else
        %In this case, the constraint is not ok, xk is not feasible. we set
        %the step direction as the derivative of the constraint
        gc = grad_const(xk);
        g = gc(:,ind);
        
        %stepsize = abs(cmax)*2;
        stepsize = abs(cmax)/norm(g)^2;
        
        xk = xk-stepsize*g;
        
        
        K = [K -1];
        
    end
    
    
end

x = get_x_real(xk);

subplot(2,1,1)
plot(x(1:128))
subplot(2,1,2)
plot(x_sol)

%figure
%plot(J)