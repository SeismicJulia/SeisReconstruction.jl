
function CGLS(x0,y,operators, parameters, μ, Ni, tol)
        
    x = copy(x0);
    r= y - LinearOperator(x,operators,parameters,adj=false);
    s =  LinearOperator(r,operators,parameters,adj=true)  - μ*x;
    p=copy(s);
    J=zeros(Ni)

    gamma= InnerProduct(s,s);
    norms0=sqrt(gamma); #norm of the gradient is used to stop.
    k=0;
    flag=false;
    
    
    while k <Ni && flag == false ;


        k = k+1;

        println("Nested iteration i= $k")
        
        q = LinearOperator(p,operators,parameters,adj=false);
        delta =  InnerProduct(q,q) + μ*InnerProduct(p,p);
        if delta == 0; #to avoid divistion by zero
             delta = 1.e-10;
        end

        alpha = gamma/delta; 
        x = x + alpha*p; #update the model
        r = r - alpha*q;
        s = LinearOperator(r,operators,parameters,adj=true) - μ*x;
        gamma1  = InnerProduct(s,s);
        norms  = sqrt(gamma1);
        beta = gamma1/gamma;
        gamma = gamma1;
        p = s + beta*p;
        flag = (norms<=norms0 * tol);
        
        if flag == true
            println("CGLS convergence reached at k= $k iterations.")
        end
        
        
        nres = norms / norms0;
        
        e = LinearOperator(x,operators,parameters,adj=false)-y; 
        J[k] = (sum((abs.(e[:]))).^2) + μ*sum((abs.(x[:])).^2 );

        if k == Ni
            println("CGLS reached the maximum number of iterations Nmax= $Ni.")
        end


    end


    return x, J

end
