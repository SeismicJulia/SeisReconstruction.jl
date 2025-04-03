using SeisPlot
using SeisReconstruction
using SeisProcessing
using LinearAlgebra
using DSP
using PyPlot

dt=0.002; # dt 
nt=512; #number of simples in time dimension
nx1=64; # number of samples in x1 dimension
nx2=64; # number of samples in x2 dimension
dx1=10.0; 
dx2=10.0;
f0= 25.0; # dominant frequency for the wavelet
tau=[0.256, 0.512, 0.768]; # t0 for the events
amp=[1.0, 0.5, -1.0]; # amplitude 
p2=[0.3,0.1, 0.45] ; # ray parameter dimension 1
p1=[0.15,0.2, 0.15]; # ray paramteter dimension2 



dtrue=SeisParabEvents( tau=tau,amp=amp, p2=p2,p1=p1, dt=dt, nt=nt, dx1=dx1, nx1=nx1, dx2=dx2, nx2=nx2, f0=f0)*(10^2);
dtruen=SeisAddNoise(dtrue, 1.0, db=false, L=6);


nt,nr,ns=size(dtruen);

dobsn =copy(dtruen); # scale

#dobsn= SeisAddNoise(dobs, .8, db=false, L=6);


dobsn=SeisDecimate(dobsn;mode="random",perc=40);
S = CalculateSampling(dobsn)

patch_size=(32,16,16); #Patch size in LocalFourier Operator5mu=μ,to
Noverlap=(16,8,8); #Overlap of patches in LocalFourier Operator
dims=size(dobsn); #Parameter to recover the right dimensions.
μ= 5.0;
tolout=1.e-10;
parameters=  [Dict(:w => S),
Dict(:patch_size=>patch_size, :Noverlap=>Noverlap, :dims=>dims, :normalize=>true, :padd=>false)];
operators=[WeightingOp, LocalFFTOp];
x0 = zeros(ComplexF64,spec_size(dobsn,patch_size,Noverlap));
m, J = FISTA(x0,dobsn,operators,parameters, μ=μ,Ni=200,tol=tolout);



drec=real(LocalFFTOp(m,false; patch_size, Noverlap, dims, normalize=true, padd=false));


difference = drec .- dtrue;

# Example slice index and data (replace with your actual data)
ix = 32  # Inline or crossline index

# Create figure
fig = PyPlot.figure(fignum=1, figsize=(15, 12))
# Panel 1: dtrue
subplot(141)
SeisPlotTX(dtrue[:, :, ix], style="wiggles", pclip=95, dx=25.0, dy=0.004, fignum=1, wiggle_trace_increment=3)
title("dtrue")
xlabel("Offset (m)")
ylabel("Time (s)")
gca().tick_params(labelsize=10)

# Panel 2: dobs
subplot(142)
SeisPlotTX(dobsn[:, :, ix],style="wiggles", pclip=95, dx=25.0, dy=0.004, fignum=1, wiggle_trace_increment=2)
title("dobs")
xlabel("Offset (m)")
ylabel("Time (s)")
gca().tick_params(labelsize=10)

# Panel 3: drec
subplot(143)
SeisPlotTX(drec[:, :, ix], style="wiggles", pclip=95, dx=25.0, dy=0.004, fignum=1, wiggle_trace_increment=2)
title("drec")
xlabel("Offset (m)")
ylabel("Time (s)")
gca().tick_params(labelsize=10)

# Panel 4: difference
subplot(144)
SeisPlotTX(difference[:,:,ix], style="wiggles", pclip=95, dx=25.0, dy=0.004, fignum=1, wiggle_trace_increment=2)
title("drec")
title("difference")
xlabel("Offset (m)")
ylabel("Time (s)")
gca().tick_params(labelsize=10)

tight_layout()
