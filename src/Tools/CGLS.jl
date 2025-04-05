
function CGLS(m0,dobs,operators, parameters; μ::Union{AbstractFloat,Vector}=0.5, Ni=100, tol=1e-6, history=true)


    if history


        println("")
        println(" =====================================================================================")
        println("                      Conjugate Gradients Least Squares (CGLS)")
        println(" =====================================================================================")
        println("")
    
    end




    if μ isa AbstractFloat
        μ = fill(μ, Ni)  # Create a vector of size Ni
    elseif μ isa Vector
        μ = float.(μ)  # Ensure it's Float64
    else
        error("μ should be either a Float or a Vector.")
    end
    

    #Initialize history tracking
    J_history = Float64[]
    Jm_history = Float64[]
    Jr_history = Float64[]
    grad_norm_history = Float64[]

    #Initial models and parameters that can change
    m = m0;
    y=dobs;
    Jm0=norm(dobs,2)^2;
    Jr0=norm(m0,2)^2;
    J0 = Jm0 + μ[begin]*Jr0
    k=0; #Initialize couter


    if history
        header = "k         ||y-Ax||²₂              ||x||²₂                   μ                   J"
        println(""); 
        println(header);
      #  @printf("%3.0f %20.10e %20.10e  %20.10e %20.10e\n",0,Jm0, Jr0,μ[begin],J0);
    end

    r = y .- LinearOperator(m,operators,parameters,adj=false);  #residual in data space
    ∇J = LinearOperator(r,operators,parameters,adj=true) - μ[begin]*m #Gradient in the model space
    d = ∇J # First conjugate direction





    #Main loop

    
    while k < Ni

        k += 1  # Increment iteration counter


        Ad=LinearOperator(d,operators,parameters,adj=false);

        sl_denom = InnerProduct(Ad, Ad) + μ[k]*InnerProduct(d, d)  #Step length denom

        if sl_denom < tol
            println("sl_denom ≈ 0. It is not possible to compute the step length to move in the conjugate direction!. Loop ending at iteration ",k)
            break
        end

        α = InnerProduct(∇J, ∇J) / sl_denom;# Step length for the conjugate direction
        
        m_new = m + α * d  # Update model parameters
        r = y - LinearOperator(m_new, operators, parameters, adj=false) # Update residual
        ∇J_new=  LinearOperator(r,operators,parameters,adj=true) - μ[k]*m_new   # Update gradient        
        
        β = InnerProduct(∇J_new, ∇J_new) / InnerProduct(∇J,∇J)  # Compute β for CG update
        d = ∇J_new + β*d  # Update conjugate direction
    
        Jmk= norm(r)^2; #New misfit
        Jrk= norm(m_new)^2; #New model norm
        Jk = Jmk +μ[k]*Jrk; # New cost function value

        res_norm = norm(r, 2)
        grad_norm= norm(∇J_new,2);


        push!(Jm_history, Jmk); #save misfit at each iteration
        push!(Jr_history, Jrk); #save regularization at each iteration
        push!(J_history,Jk); #save objective at each iteration
        push!(grad_norm_history, grad_norm) # save grad norm at each iteration




        #Tolerance cheking and printing


        if history && k <Ni
            @printf("%3.0f %20.10e %20.10e  %20.10e %20.10e\n", k, Jmk, Jrk, μ[k], Jk)
        end        


        if res_norm < tol && grad_norm < tol

            if history
                println("CGLS converged at iteration $k with residual norm: $res_norm and gradient norm: $grad_norm .")
            end
            break
        end

        if length(J_history) > 1 && J_history[end] > eps()

            ΔJ= abs((J_history[end] - J_history[end-1]) / J_history[end])

            if round(ΔJ,digits=8) < tol
               println("Loop for CGLS stopped at $k iterations.")
               println("REASON: ")
               println(" ΔJ = $ΔJ   is < than the established tolerance = $tol used.")
               break
            end			
        end

        #Update model and gradient for residual and step length iterations

        m= m_new;
        ∇J= ∇J_new;


    end

        #possible outputs

        #misfit_iteration= res_norm; # residual norm or misfit => ensure you fit the data
        #grad_norm_history # norm of ∇J => ensure you reached a stable solution
        #Jm_history; # data fidelity term for χ² test.
        #Jr_history; # model norm fidelity term for χ² test.
 

    return m, J_history




end    