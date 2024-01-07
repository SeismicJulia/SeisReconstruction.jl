function FISTA(x0,y,operators,parameters; μ= 0.5,Ni=100,tolerance=1.0e-3,history=true)

    println("")
    println(" ==================================================")
    println(" Fast Iterative Soft Thresholding Algorithm (FISTA)")
    println(" ==================================================")
    println("")

    #Initialize
    x0=randn(size(x0));
    α = PowerMethod(x0,operators,parameters);
    η= 0.95/α;
    m = zeros(Float64,size(x0));
    t=1.0;
    p = copy(m);
    J0=norm(y,2)^2;
    misfit_term= 0.0;
    regularization_term= 0.0;
    norm_cost=1.0;
    Je=Float64[]
    push!(Je,norm_cost)
    k=0;


    
    if history
        header = " k          ||y-Ax||²₂                ||x||₁               μ              J/J0"
        println(""); 
        println(header);
        @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k,misfit_term, regularization_term,μ, norm_cost);

     end
    
    while k < Ni

            #Update counter
            k=k+1;

            #move and gradient direction and threshold
            m_old = m; # save model to perturb;
            Am=LinearOperator(p,operators,parameters,adj=false);
            r= Am .-y; #residual= (A*x.-y)
            ∇fₖ=LinearOperator(r,operators,parameters,adj=true); #A'(r)
            m = p .- η*∇fₖ; #update in GD direction.
            m = SoftThresholding.(m,η,μ); #soft-thershold the  update
            
            ##FISTA acceleration step
            t_old=t;
            t = (0.5)*(1.0 + sqrt(1.0 + 4.0*(t_old)^2));
            p = m +((t_old-1.0)/t)*(m-m_old);

            yp= LinearOperator(m,operators,parameters,adj=false);
            misfit_term= sum(abs.(yp .- y).^2);
            regularization_term= sum(abs.(m));  
            cost = (1/2)*misfit_term + μ*regularization_term; 
            norm_cost = cost/J0; 
            push!(Je,norm_cost);
            
            if history
                @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k, misfit_term, regularization_term,μ, norm_cost);
            end
            
            if k> 1
                ΔJ= abs(Je[k] - Je[k-1])/((Je[k]+Je[k-1])/2);
                if ΔJ < tolerance
                println("Loop ended at $k iterations.")
                println("REASON: ")
                println(" ΔJ = $ΔJ  is < than the established tolerance = $tolerance used.")
                break
                end
            end

    end
        return m, Je
end
