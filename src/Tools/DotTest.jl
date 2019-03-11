"""
    DotTest(m_rand,d_rand,operators,parameters)
Dot product test for a vector of linear operators
See also: [`ConjugateGradients`](@ref)

"""
function DotTest(m_rand,d_rand,operators,parameters)

	m_adj = LinearOperator(d_rand,operators,parameters,adj=true)
	d_fwd = LinearOperator(m_rand,operators,parameters,adj=false)
	inner1 = InnerProduct(d_rand,d_fwd)
	inner2 = InnerProduct(m_rand,m_adj)
	println("<d_rand,d_fwd> = ",inner1)
	println("<m_rand,m_adj> = ",inner2)
	return inner1,inner2

end
