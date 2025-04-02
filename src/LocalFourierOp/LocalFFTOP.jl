

function LocalFFTOp(in, adj::Bool; patch_size, Noverlap, dims,normalize, padd)
    
    norm = normalize ? sqrt(length(in[:])) : 1.0 # Normalize
    
    #adj= true compute the short time fourier transform of the signal
    if adj==true
        
        out= stfft(in,patch_size,Noverlap,padd);   
        return out
      
    else

        nt,nx,ny=dims;
        out, θ= istfft(in,patch_size, Noverlap, padd; Nt=nt, Nx=nx, Ny=ny);
        out= (out.*θ);
        return out#, θ #θ It should give also theta, or should find a way of computing theta outside?
    end
end
