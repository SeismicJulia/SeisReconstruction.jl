"""
    ConjugateGradients(d,operators,parameters;<keyword arguments>)
Conjugate Gradients following Algorithm 2 from Scales, 1987.
The user provides an array of linear operators. Verify that linear operator(s) pass the dot product.
See also: [`DotTest`](@ref)

# Arguments
- `Niter=10` : Number of iterations
- `mu=0`
- `tol=1.0e-15`

"""
function ConjugateGradients(d,operators,parameters;Niter=10,mu=0,tol=1.0e-15)

    cost = Float64[]
    r = copy(d)
    g = LinearOperator(r,operators,parameters,adj=true)
    m = zero(g)
    s = copy(g)
    gamma = InnerProduct(g,g)
    gamma00 = gamma
    cost0 = InnerProduct(r,r)
    push!(cost,1.0)
    for iter = 1 : Niter
	t = LinearOperator(s,operators,parameters,adj=false)
	delta = InnerProduct(t,t) + mu*InnerProduct(s,s)
	if delta <= tol
#	    println("delta reached tolerance, ending at iteration ",iter)
	    break;
	end
	alpha = gamma/delta
	m = m + alpha*s
	r = r - alpha*t
	g = LinearOperator(r,operators,parameters,adj=true)
	g = g - mu*m
	gamma0 = copy(gamma)
	gamma = InnerProduct(g,g)
        cost1 = InnerProduct(r,r) + mu*InnerProduct(m,m)
        push!(cost,cost1/cost0)
	beta = gamma/gamma0
	s = beta*s + g
	if (sqrt(gamma) <= sqrt(gamma00) * tol)
	    println("tolerance reached, ending at iteration ",iter)
	    break;
	end
    end

    return m, cost
end
