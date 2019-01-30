using SeisProcessing
using Base.Test

# test that FFTOp passes the dot product test

m_rand = rand(ComplexF32,(20,20,20,20));
d_rand = rand(ComplexF32,(20,20,20,20));

T = rand(20,20,20,20);

a = (LinearIndices(T))[findall(T.<=0.5)];
b = (LinearIndices(T))[findall(T.>0.5)];

T[a] .= 0.;
T[b] .= 1.;


operators = [WeightingOp,FFTOp,WeightingOp]
parameters = [Dict(:w=>T),Dict(:normalize=>true),Dict(:w=>fill(1,size(T)))]

a,b = SeisReconstruction.DotTest(m_rand,d_rand,operators,parameters)
@test abs((a-b)/a) < 0.01

