@testset "FFTOp" begin
    # test that a linear operator passes the dot product test.
    # See, for instance, Earth Soundings Analysis: Processing Versus Inversion by Jon Clearbout

    # The lineaer operaor is given by L =  T FFTop W
    # where T is sampling,FFTop is an inverse FFT and W is applying weights

    m_rand = rand(ComplexF64,(20,20,20,20));
    d_rand = rand(ComplexF64,(20,20,20,20));

    T = rand(20,20,20,20);

    a = (LinearIndices(T))[findall(T.<=0.5)];
    b = (LinearIndices(T))[findall(T.>0.5)];

    T[a] .= 0.;
    T[b] .= 1.;


    operators = [WeightingOp,FFTOp,WeightingOp]
    parameters = [Dict(:w=>T),Dict(:normalize=>true),Dict(:w=>fill(1,size(T)))]

    a,b = DotTest(m_rand,d_rand,operators,parameters)
    @test abs((a-b)/a) < 0.01
end
