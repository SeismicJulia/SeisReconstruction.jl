
function PowerMethod(x0, operators, parameters; tol=1e-2, Ni=50)




    x = normalize(x0)  # Normalize initial vector
    λ = 0.0 #initialize eigen value

    k=0;

    while k  < Ni


        k=k+1

        println("Power Iteration: $k")

        aux = LinearOperator(x, operators, parameters, adj=false)  # Applying A
        y = LinearOperator(aux, operators, parameters, adj=true)  # Applying A^H A

        λnew = norm(y, 2)  # Compute new eigenvalue
        xnew = y / λnew  # Normalize new eigenvector

        # Convergence check
        if abs(λnew - λ) < tol
            println("Power iteration converged after $k iterations.")
            return λnew, xnew
        end

        λ =λnew
        x = xnew

    end

    println("Power iteration reached max iterations ($Ni) without full convergence.")
    
    return λ, x
end



function MinimumEigenValue(x0, λmax, vmax, operators, parameters; num_tests=10)
    

    println("Estimating λmin from λmax and vmax.")

    # Compute the dominant eigenvector using Power Iteration
    
    #λ_max, v_max = PowerMethod(x0, operators, parameters; tol=tol, max_iter=max_iter)

    λmin = Inf  # Initialize with a large value

    for k in 1:num_tests;


        println(" Iteration = $k")
        v = randn(size(x0)) # Generate a random vector to test against vmax
        v = v .- (InnerProduct(v, vmax)/InnerProduct(vmax, vmax)).*vmax;  #Make it orthogonal to vmax
        v = normalize(v)  # Normalize
        Av = LinearOperator(v, operators, parameters, adj=false)
        AᴴAv = LinearOperator(Av, operators, parameters, adj=true)  # Compute A^H A v
        λ = InnerProduct(v,AᴴAv)/InnerProduct(v,v) #estimate the correspondent minimum eigen value
        λmin = min(λmin, λ)
    end

    return λmin
end





#=
function PowerMethod(x0,operators,parameters)
    
    x= x0;
    α=0.0;
    
    begin 
            for k = 1:20;
            
            println("Iteration: $k")
            
            aux=LinearOperator(x,operators,parameters,adj=false)
            y=LinearOperator(aux,operators,parameters,adj=true)
            n = norm(y,2);
            x = y/n;
         α = n;
        end
    end

    return α
end





function PowerMethodWithRayleigh(x0, operators, parameters)
    x = x0
    α = 0.0

    for k = 1:5
        println("Power iteration number $k")
        
        # Apply the linear operator (adj=false) to get the auxiliary vector
        aux = LinearOperator(x, operators, parameters, adj=true)
        
        # Apply the adjoint of the linear operator (adj=true)
        y = LinearOperator(aux, operators, parameters, adj=false)
        
        # Normalize the transformed vector
        n = norm(y, 2)
        x = y / n
        
        # Compute the Rayleigh quotient: λ = (x' * A * x) / (x' * x)
        # Here, 'A' is represented by the linear operator applied twice as above.
        num = dot(x, LinearOperator(x, operators, parameters, adj=false))  # x' * A * x
        denom = dot(x, x)  # x' * x
        
        α = num / denom  # Rayleigh quotient as the eigenvalue approximation

        println("Rayleigh quotient eigenvalue estimate: $α")
    end

    return α
end






#=

function PowerMethodWithRayleighNonSquare(x0, operators, parameters)
    
    
    x = x0
    α = 0.0

    for k = 1:20
        
        println("Power iteration number $k")

        # Apply the linear operator A (adj=false) to x
        y = LinearOperator(x, operators, parameters, adj=false)
        
        # Apply the adjoint of the linear operator A^* (adj=true) to y
        z = LinearOperator(y, operators, parameters, adj=true)
        
        # Normalize the transformed vector z
        n = norm(z, 2)
        x = z / n  # Normalizing to avoid overflow
        
        # Here we cannot compute the Rayleigh quotient as (x' A x), 
        # because A maps between different spaces.
        # Instead, we use z' * A * x (since z is in the appropriate space after adjoint application)
        num = dot(z, LinearOperator(x, operators, parameters, adj=false))  # z' * A * x
        denom = dot(z, z)  # z' * z (to stay in the right space dimension)
        
        α = num / denom  # Generalized Rayleigh quotient

        println("Generalized Rayleigh quotient eigenvalue estimate: $α")
    end

    return α
end






function PowerMethodWithRayleigh(x0, operators, parameters)
    x = x0
    α = 0.0

    for k = 1:20
        println("Power iteration number $k")
        
        # Apply the linear operator (adj=false) to get the auxiliary vector
        aux = LinearOperator(x, operators, parameters, adj=true)
        
        # Apply the adjoint of the linear operator (adj=true)
        y = LinearOperator(aux, operators, parameters, adj=false)
        
        # Normalize the transformed vector
        n = norm(y, 2)
        x = y / n-
        
        # Compute the Rayleigh quotient: λ = (x' * A * x) / (x' * x)
        # Here, 'A' is represented by the linear operator applied twice as above.
        num = dot(x, LinearOperator(x, operators, parameters, adj=false))  # x' * A * x
        denom = dot(x, x)  # x' * x
        
        α = num / denom  # Rayleigh quotient as the eigenvalue approximation

        println("Rayleigh quotient eigenvalue estimate: $α")
    end

    return α
end

function power_iteration_non_square(A; tol=1e-6, max_iters=1000)
    m, n = size(A)  # Get the size of the matrix (non-square allowed)

    # Start with a random vector of appropriate size
    if m > n
        b = rand(ComplexF64, n)  # Start with a vector in the domain size (n)
    else
        b = rand(ComplexF64, m)  # Start with a vector in the codomain size (m)
    end
    
    # Normalize the starting vector
    b = b / norm(b)
    
    λ_old = 0.0 + 0.0im  # Initialize the eigenvalue guess as a complex number
    for i in 1:max_iters
        # Apply the matrix to the vector
        b_new = A * b
        
        # Apply the adjoint (conjugate transpose) to handle non-Hermitian matrices
        z_new = A' * b_new  # Conjugate transpose (adjoint operator)
        
        # Normalize the resulting vector
        b_new = b_new / norm(b_new)
        
        # Estimate the eigenvalue using the Rayleigh quotient (z_new' * A * b_new)
        λ = (b_new' * A * b_new) / (b_new' * b_new)  # λ ≈ Rayleigh quotient for complex matrices
        
        # Check for convergence (absolute difference of eigenvalue estimates)
        if abs(λ - λ_old) < tol
            println("Converged after $i iterations")
            return λ, b_new
        end
        
        # Update the vector and eigenvalue for the next iteration
        b = b_new
        λ_old = λ
    end
    
    println("Maximum iterations reached")
    return λ_old, b
end
=#


=#