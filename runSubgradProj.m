function [xsolved] = runSubgradProj(MAXIT,step)

    load('cs.mat');

    %declarations
    x_sol = x;
    Fopt=0;
    it = 0;

    %Decomposition of the problem in real and imaginary parts + null imaginary
    %part constraint
    X_us2 = [real(X_us); imag(X_us); zeros(128,1)];
    F_us2 = [real(F_us) -imag(F_us); imag(F_us) real(F_us); zeros(128,128) eye(128)];

    %--------------------------------------
    % min norm(x,1) 
    % s.t. F_us*x - X_us = 0
    %      imag(x) = 0
    %      real(x) >= 0
    %
    %--------------------------------------

    %Objective fuction
    func = @(x) norm(x(1:128) + 1i*x(129:256), 1);

    % Equality constraint set projector
    PROJ = pinv(F_us2);

    %sub gradient of norm(x,1) = sign(x).
    grad = @(z) sign(z);

    % Constraint (Ax - b <= 0)
    const = @(z) -[z(1:128); zeros(128,1)];
    grad_const = @(z) -[ones(128,1); zeros(128,1)];

    f = [];

    %Starting point : LS solution
    xk = pinv(F_us2)*X_us2;

    while(it < MAXIT)
        it = it+1;

        %Finding most violated contraint, should be positive if violated
        [cmax, ind] = max(const(xk));

        if(cmax <= 0)
            %In this case, the inequality constraint is not violated.
            %We use the normal proj subgrad method

            g = grad(xk);
            g = g - ([zeros(128,1); 10*ones(128,1)] == g)/10; %Putting zeros to -0.1 to enforce positivity

            %stepsize
            stepsize=step/it;

            %Update
            xk = xk - stepsize*(eye(256) - PROJ*F_us2)*g; 

        else
            %In this case, the constraint is not ok, xk is not feasible. we set
            %the step direction as the derivative of the constraint
            g = zeros(256,1);
            gc = grad_const(xk);
            g(ind) = gc(ind);

            stepsize = abs(cmax)/norm(g)^2;

            xk = xk-stepsize*g;

        end

        %Keeping objective function value
        objective(it)=func(xk)-Fopt;

    end

    xsolved=xk;

    %Plotting convergence graph
    semilogy((objective(2:end)),'LineWidth',1.5)
    xlabel('Number of Iterations');ylabel('f(x_k)');
    grid on
end
