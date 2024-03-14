function ADMM(m0,d_obs,operators,parameters; ρ= 0.5, μ= 5.0 ,tolin=1e-2, tolout=1e-4, Ni=5,Ne=50, history=true)
    
    
    #Initialize variables and allocate arrays:

    J0=norm(d_obs,2)^2; #Initial objective function value
    misfit_term=0.0; regularization_term= 0.0 ;
    norm_cost=1.0;
    Je=Float64[]; #External cost function 

    u=zeros(size(m0)); 
    z=zeros(size(m0)); 
    w=zeros(size(m0)); 
    
    push!( Je,norm_cost)

    k=0; # Outer loop counter

     
    if history
        header = " k          ||y-Ax||²₂                ||x||₁               μ              J/J0"
        println(""); 
        println(header);
        @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k,misfit_term, regularization_term,μ, norm_cost);
      
     end

    while k < Ne;  

    
        #ADMM
        k=k+1 #Update counter    
        b=  -1*LinearOperator(z.- u ,operators, parameters, adj=false) .+ d_obs; # Change variables
        #wcgls, Jcgls = ConjugateGradients(b,operators,parameters;mu=ρ,tol=tolin, Niter=Ni);
        w, Jcgls =CGLS(m0,b,operators, parameters, ρ, Ni, tolin); #CGLS ==> x-update
        x=  w .+ (z .- u) ; #change variables
        z= SoftThresholding.( x.+ u,ρ,μ) ; #Soft-thresholding ==> z-update
        u= u .+ (x .-z); #Dual update

        #Objective function

        yp= LinearOperator(z,operators,parameters,adj=false); #predicte
        misfit_term= sum(abs.(yp .- d_obs).^2);  #|| Ax - y ||₂²
        regularization_term= sum(abs.(z));  #|| x ||₁
        cost = (1/2)*misfit_term + μ*regularization_term; #loss function
        norm_cost = cost/J0; #normalize costfunction
        push!(Je,norm_cost);
    
        #Print information:
        if history
            @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k, misfit_term, regularization_term,μ, norm_cost);
         end

        #Tolerance
        if k> 1
            ΔJ= abs(Je[k] - Je[k-1])
            if ΔJ < tolout
               println("Loop ended at $k iterations.")
               println("REASON: ")
               println(" ΔJ = $ΔJ  is <than the established outer loop tolerance = $tolout used.")
               break
            end			
        end
    end

    println("Outer-loop ended at $k iterations.")

    
    return z, Je

end