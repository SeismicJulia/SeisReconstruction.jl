#Careful! needs LinearAlgebra in Reconstruction.jl
"""
 SeisSpitzFX(d;<keyword arguments>)

 First order f-x interpolation using Spitz' method. 
 Input size is size(d)= (nt,nh). Output size is (nt,2*nh-1)

# Arguments
-`dt`: sampling interval in secs
-`npf`: prediction filter length
- `pre1` : percentage of pre-whitening for Pef estimation
- `pre2` : percentage of pre-whitening for estimation of missing samples
- `flow` : minimun temporal frequency to reconstruct (Hz)
- `fhigh` : maximun temporal frequency to reconstruct (Hz)

# Example
'''julia
d_interp = SeisSpitzFX(d)
'''
*  Reference: Spitz, S., 1991, Seismic trace interpolation in the F-X domain: Geophysics, 56, 785-794.
"""
function SeisSpitzFX(d;dt=0.004,npf=2,pre1=1,pre2=1,flow=1,fhigh=90)


# 1) Fourier Transform to FX domain 

 vect = size(d);
 nt = vect[1]
 nf1 = nextpow(2,vect[1]);
 nf2 = 2*nf1 
 nh = vect[2]


 d1,ifl1,ifh1 = pad(d,dt,flow,fhigh)
 D1 = fft(d1,1)

 d2,ifl2,ifh2 = pad(d,dt,flow,fhigh;npad=2)
 D2 = fft(d2,1)


 INTDF = zeros(typeof(D2[1]),nf1, 2*nh-1);

 for ia = ifl1:ifh1;

# 2) Select frequency slices from frequencies f and f/2

    x1 = transpose(D1[ia,:]);
    x2 = transpose(D2[ia,:]);
# 3) Estime the PFs needed for frequency sample ia;

    PF = prediction_filter(x2,npf,pre1);

# 4) Interpolate spatial data at freq. sample ia
    y = interpolate_freq(x1,PF,pre2);

# 5) Replace interpolated frequency slice into final matrix

    INTDF[ia,:] = transpose(y);

end

# 6) Honor symmetries

nf1_2 = convert(Int,round(nf1/2))

INTDF[nf1_2+2:nf1,:]=conj(reverse(INTDF[2:nf1_2,:],dims=1));

# 7) Antitransform from f-x to t-x
 d_interp = real(ifft(INTDF,1));

 d_interp = d_interp[1:nt,:];

 return d_interp
end

#*******************************************************************************
"""
**INTERPOLATE_FREQ**
*fx_interpolation for reconstruction of traces*

**IN**
*   x:   spatial data to be recostruction at a given freq.
*       PF:  prediction error filter
*       pre: prewhitening in percentage

**OUT**
*  y:   reconstructed vector of data at a particular frequency
"""
function interpolate_freq(x,PF,pre)


 np = length(PF);
 nx = length(x);
 ny = 2*nx-1;       # length of interpolated signal

# Forward step using PFs for interpolation
 TMPF1 = hcat(reverse(transpose(PF),dims=2),-1)
 W1 = transpose(Convmtx(TMPF1,ny-np));


# Backward step

 TMPF2 = conj(reverse(TMPF1,dims=2));
 W2 = transpose(Convmtx(TMPF2,ny-np));

# Forward and backward combined in an augmented system

 WT = [W1;W2]

# Separation of  Matrix into  known and unknown parts
 A = WT[:,2:2:ny];

B = -1 .* WT[:,1:2:ny]*transpose(x);
# Least squares solution for missing data, prewhitening
# is added to guarantee stability

 R= A'*A;
 g = A'*B;


 mu = (pre/100.)*tr(R)/(nx-1);

y1 =  (R+mu*Matrix(I,nx-1,nx-1))\g;

 y = zeros(typeof(x[1]),1,ny);

y[1:2:ny]=x;
 y[2:2:ny]=transpose(y1);

 return y
 end

 #******************************************************************************
"""
**PREDICTION_FILTER**
* compute prediction filters

**IN**
*   VEC:   spatial data at a given frequency
*       np:    length of pef
*       pre:   prewhitening in percentage

**OUT**
*  PF:    prediction filter
"""
function prediction_filter(xi,np,pre);

 ns = length(xi);
 

C = zeros(typeof(xi[1]),ns-np+1,np+1);

 for j=1:ns-np

    C[j,:]=xi[j+np:-1:j];

 end
 
 # Make kernel for pef estimation (forward and backward prediction included)

 A = [C[:,2:np+1];conj(reverse(C[:,1:np],dims=2))];

# RHS of system of equations

 B = [C[:,1];conj(C[:,np+1])];

 # Solve A.PF = B using least-squares with pre-whitening

 R = A'*A;
 g = A'*B;

 mu = real((pre/100.)*tr(R)/np);

 PF =  (R+mu*Matrix(I,np,np))\g;   # same as inv(R)*g

 return PF
 end


#*****************************************************

function pad(d,dt,flow,fhigh;npad=1)


vect = size(d)
nt = vect[1]

#calculate index frequencies

nf = npad*nextpow(2,nt)

dpad = zeros(nf,vect[2])

dpad[1:nt,:] = d 

ilow = convert(Int64,floor(flow*nf*dt)+1)
ihigh = convert(Int64,floor(fhigh*nf*dt)+1)

return dpad,ilow,ihigh
end
