function [x, it, prec] = backtracking_gradient_descent(x0, alpha, beta, f, grad)

    x = x0;
    it = 0;
    
    tol = e-3;
    
    grad_x = feval(grad, x);
    prec = norm(grad_x);
    
    while(it < 100 || prec > tol)
        
        
       
        %Step direction
        delta_x = -grad_x;
         
        
        %Backtracking line search
        t = 1;
        while(feval(f, x+t*delta_x) > feval(f, x) + alpha*t*grad_x'*delta_x)
           
            t = beta*t;
            
        end
        
        
        %Update
        x = x + t*delta_x;
        grad_x = feval(grad, x);
        prec = norm(grad_x);
        
        it = it+1;
        
    end

end