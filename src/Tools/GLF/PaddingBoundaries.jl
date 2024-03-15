

function padd_boundaries(in,Nw, Nkx, Nky)
    
    paddt = zeros(Float64,(Nw ÷ 2))
    paddx = zeros(Float64,(Nkx ÷ 2)) # Pad to extent function to have an integr number of window and include all the samples.
    paddy = zeros(Float64,(Nky ÷ 2))

    out1=zeros(typeof(in[1,1,1]),(length(paddt)+size(in,1) +length(paddt)), (length(paddx)+size(in,2) +length(paddx)), (length(paddy) + size(in,3) + length(paddy)))

    out1[length(paddt)+1 : size(out1,1)-length(paddt),length(paddx)+1: size(out1,2)-length(paddx),length(paddy)+1: size(out1,3)-length(paddy)]= in[:,:,:]
  
    
    
    return out1;
    
end



function padd_recover_amp(in,Nw,Nkx,Nky,hopt,hopx,hopy)
    
    Nt= size(in,1); Nx=size(in,2); Ny=size(in,3);
    paddt= zeros(Float64,mod(mod(-(Nt-Nw),hopt),Nw)) #Padd to recover the amplitud
    paddx= zeros(Float64,mod(mod(-(Nx-Nkx),hopx),Nkx))
    paddy= zeros(Float64,mod(mod(-(Ny-Nky),hopy),Nky))

    out2= zeros(typeof(in[1,1,1]),(Nt+length(paddt),Nx +length(paddx), Ny + length(paddy)))

    out2[1:Nt,1:Nx,1:Ny] .= in;
    Nt2,Nx2,Ny2= size(out2)

    return out2, size(out2,1), size(out2,2), size(out2,3)
end