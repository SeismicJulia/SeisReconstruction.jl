

function FISTA(x0,y,operators,parameters; μ::Union{AbstractFloat,Vector}=0.5,Ni=100,tol=1.0e-3,history=true)

    println("")
    println(" ===============================================================")
    println("      Fast Iterative Soft Thresholding Algorithm (FISTA)       ")
    println(" ==============================================================")
    println("")
    
    
    if typeof(μ) <: AbstractFloat;
        μ=μ*ones(Ni);
    elseif typeof(μ) <: Vector;
       μ=convert.(Float64,μ);
    else
       error("μ should be either AbstractFloat or Vector type.")
    end

    #Initialize
    x0=randn(eltype(x0),size(x0));
    α, vmax = PowerMethod(x0,operators,parameters);
    η= 0.95/α;
    m = similar(x0);
    t=1.0;
    p = copy(m);
    #J0=norm(y,2)^2;
    
    
 
    J_history=Float64[]
    Jr_history=Float64[]
    Jm_history=Float64[]
    grad_norm_history = Float64[]
    res_norm_history = Float64[]
    get_time=Float64[]





     
    if history
        header = " k          ||y-Ax||²₂                ||x||₁               μ              J"
        println(""); 
        println(header);
     end


     #Main Loop:

    k=0; #iteration counter
    while k <  Ni

         gt=@elapsed begin
            #Update counter
            k=k+1;
            
          #  println("-----------------------------------------------")
           # print("Monitoring iteration number and cost function "); @show k; @show J[k];
       
            #move and gradient direction and threshold
            m_old = m; # save model to perturb;
            Am=LinearOperator(p,operators,parameters,adj=false);
            r= Am .-y; #residual= (A*x.-y)
            ∇J=LinearOperator(r,operators,parameters,adj=true); #A'(r)
            m = p .- η*∇J; #update in GD direction.
            m = SoftThresholding.(m,η,μ[k]); #soft-thershold the  update
            
            ##FISTA acceleration step
            t_old=t;
            t = (0.5)*(1.0 + sqrt(1.0 + 4.0*(t_old)^2));
            p = m +((t_old-1.0)/t)*(m-m_old);
        end
            
        
        
        
            ypred= LinearOperator(m,operators,parameters,adj=false);
            Jmk= norm(y .-ypred,2)^2;#sum(abs.(yp .- y).^2);
            Jrk= sum(abs.(m));  
            Jk = Jmk + μ[k]*Jrk; 
            res_norm= sqrt(Jmk);
            grad_norm= norm(∇J,2);
            
            
            push!(J_history,Jk);
            push!(Jm_history,Jmk);
            push!(Jr_history,Jrk);
            push!(grad_norm_history,grad_norm);
            push!(res_norm_history,res_norm);
            push!(get_time,gt);






            if history
                @printf("%3.0f %20.10e %20.10e %20.10e %20.10e\n", k, Jmk, Jrk,μ[k], Jk);
            end

            if res_norm < tol && grad_norm < tol

                if history
                    println("FISTA converged at iteration $k with residual norm: $res_norm and gradient norm: $grad_norm .")
                end
                break
            end    

            
            if length(J_history) > 1 && J_history[end] > eps()
                
                ΔJ= abs((J_history[end] - J_history[end-1]) / J_history[end])
                if round(ΔJ,digits=8) < tol
                    println("Loop for FISTA stopped at $k iterations.")
                    println("REASON: ")
                    println(" ΔJ = $ΔJ   is < than the established tolerance = $tol used.")
                    break
                end
            end

    end


        #Possible outputs#
    #J_history,Jk;
    #Jm_history,Jmk;
    #Jr_history;
    #grad_norm_history;
    #res_norm_history;
        
    
    return m, J_history, get_time
end
