function [numit,x] = runSubgradReduced(MAXIT,step)
    load('cs.mat');

    %declarations
    x_sol = x;
    Fopt=3;
    idx=0;
    iter=[]; %iteration of the best solution
    it = 0;

    %Decomposition of the problem in real and imaginary parts
    X_us2 = [real(X_us); imag(X_us); zeros(128,1)];
    F_us2 = [real(F_us) -imag(F_us); imag(F_us) real(F_us); zeros(128,128) eye(128)];

    %Particular solution
    xp = pinv(F_us2)*X_us2;

    %F change of basis
    F = null(F_us2);
    
    %dimension of reduced variable
    sz = size(F);
    dimz = sz(2);
    

    %--------------------------------------
    % min norm(x,1) 
    % s.t. F_us*x - X_us = 0
    %      imag(x) = 0
    %      real(x) >= 0
    % using reduced formulation x = Fz + xp
    %--------------------------------------
    
    %Real part and imaginary part masks
    rmask = [eye(128) zeros(128,128); zeros(128,128) zeros(128,128)];
    imask = [zeros(128,128) zeros(128,128); zeros(128,128) eye(128)];

    %Get real and imaginary parts
    get_x_real = @(z) rmask*(F*z + xp);
    get_x_im = @(z) imask*(F*z + xp);

    %Objective function
    func = @(z) norm(get_x_real(z)+ 1i*get_x_im(z), 1);

    %sub gradient of objective
    grad = @(z) F'*sign(F*z + xp);

    % Constraint (Ax - b <= 0)
    const = @(z) -get_x_real(z);
    grad_const = @(z) -(rmask*F)';
    
    xk = zeros(dimz,1);

    BEST_XK = xk;

    while(it < MAXIT)
        it = it+1;

        %Finding most violated contraint, should be positive if violated
        [cmax, ind] = max(const(xk));
      
        j = 0;
        
        %Making the point feasible
        while(j < 1000 && cmax >= 1e-3)
            gc = grad_const(xk);
            g = gc(:,ind);

            stepsize = abs(cmax)/norm(g)^2;

            xk = xk-stepsize*g;

            [cmax, ind] = max(const(xk));
            j = j+1;
        end

        %Now that the point is feasible, do the subgradient step 
        g = grad(xk);

        stepsize = step/it;

        %Update
        xk = xk - stepsize*g; 

        %Check for best point
        if func(xk) < func(BEST_XK)
            BEST_XK = xk;
        end

        objective(it)=func(BEST_XK)-Fopt;

        %Calculate RMSE based # of iterations given tol
        x2 = get_x_real(BEST_XK);
        x2=x2(1:128);
        RMSE(it) = sqrt(mean((x2 - x_sol).^2));

        if RMSE(it) < 1e-3
            idx=idx+1;
            iter(idx)=it;
        end

    end
    
    %output best iteration for 1e-3 error and best solution found
    numit=iter(1);
    x = get_x_real(BEST_XK);
    
    %Plotting convergence graph
    fopt=0;
    semilogy((objective(1:end)-fopt),'LineWidth',1.5)
    xlabel('Number of Iterations');ylabel('f(x_k)-f^*');
    grid on
end
 




