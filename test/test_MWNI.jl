using SeisReconstruction
using SeisProcessing
using LinearAlgebra
using Test

d = SeisLinearEvents()

#randomly decimate 50% of traces

mask = rand(1,size(d,2));
perc=50
mask[(LinearIndices(mask .< perc/100))[findall(mask .< perc/100)]] .= 0;
mask[(LinearIndices(mask .>= perc/100))[findall(mask .>= perc/100)] ] .= 1;

ddec = Array{Float64}(undef,size(d));

for it = 1 : size(d,1)
    ddec[it,:] = d[it:it,:].*mask;
end


dmwni = SeisMWNI(ddec,dt=0.004,fmax=100,Niter_external=5,Niter_internal=10)

# test that quality factor is greater than 10 Decibels
quality_factor = 10*log10(norm(d[:],2)/norm(dmwni[:]-d[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 5
