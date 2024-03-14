using SeisPlot, SeisReconstruction, SeisProcessing
using  PyPlot, LinearAlgebra, FFTW, DSP, Statistics, LinearAlgebra;


dt=0.002; # dt 
nt=512; #number of simples in time dimension
nx1=64; # number of samples in x1 dimension
nx2=64; # number of samples in x2 dimension
dx1=10.0; 
dx2=10.0;
f0= 25.0; # dominant frequency for the wavelet
tau=[0.256, 0.512, 0.768]; # t0 for the events
amp=[1.0, 0.5, -1.0]; # amplitude 
p2=[0.3,0.0, -0.45] ; # ray parameter dimension 1
p1=[0.15,0.0, -0.15]; # ray paramteter dimension2 



data=SeisParabEvents( tau=tau,amp=amp, p2=p2,p1=p1, dt=dt, nt=nt, dx1=dx1, nx1=nx1, dx2=dx2, nx2=nx2, f0=f0);

nt,nr,ns=size(data);

d_obs =(10^3)*copy(data); # scale


d_obsn= SeisAddNoise(d_obs, 1.0, db=false, L=9);


d_obsn=SeisDecimate(d_obsn;mode="random",perc=50);


S = CalculateSampling(d_obsn)

ρ=0.5;
μ= 5.0;
tol=1.e-4;
x0 = randn(Float64,size(d_obs));
operators=[ WeightingOp, FFTOp];
parameters= [Dict(:w =>S), Dict( :normalize=>true)];
m, J= ADMM(x0,d_obs,operators,parameters, ρ= ρ , μ=μ, Ne=100, Ni=100,tolin=1e-4, tolout=1e-4 );