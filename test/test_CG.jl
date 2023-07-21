@testset "CG" begin
    # test of Conjugate Gradients and compare results with closed form solution
    # Find m such that || L m - d||_2^2 + mu || m ||_2^2 is minimum

    nd = 1000
    nm = 300
    L = randn(nd,nm)
    m = randn(nm)

    # Make data and add noise 
    d = L*m
    d = d + 0.2*randn(nd)

    # closed form solution
    mu = 0.2
    m1 = (L'*L + mu*Matrix(I,nm,nm))\(L'*d)

    # now test CG

    m0 = zeros(nm)
    m2,cost = ConjugateGradients(d,[MatrixMultiplyOp],[Dict(:matrix=>L)],Niter=nm,mu=mu)

    # Compute SNR between CG solution and closed form solutions 

    quality_factor = 10*log10(norm(m1[:],2)/norm(m2[:]-m1[:],2))
    println("Quality factor = ",quality_factor)
    @test quality_factor > 10.
end
