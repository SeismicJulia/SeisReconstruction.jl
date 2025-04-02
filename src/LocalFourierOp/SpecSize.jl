

function spec_size(in,patch_size,Noverlap)
    
    Nw,Nkx,Nky= patch_size#(32,32,32)#patch_size;
    Noverlap_t, Noverlap_x, Noverlap_y = Noverlap#(16,16,16); #Number of overlap in each direction 
    hopt,hopx,hopy=check_parameters(in,Nw,Nkx,Nky,Noverlap_t,Noverlap_x, Noverlap_y; analysis_win=nothing);
    
    Nt,Nx,Ny=size(in);
    out1=padd_boundaries(in,Nw,Nkx,Nky)
    out2,Nt2,Nx2,Ny2=padd_recover_amp(out1,Nw,Nkx,Nky,hopt,hopx,hopy)

    ω = Nw;  kx=Nkx; ky=Nky 

    τ=Int(floor((Nt2-Noverlap_t)/(Nw-Noverlap_t)));
    x=Int(floor((Nx2-Noverlap_x)/(Nkx-Noverlap_x)));
    y=Int(floor((Ny2-Noverlap_y)/(Nky-Noverlap_y))); #This is for now!
    
    return (ω,kx,ky,τ,x,y)
    
end
