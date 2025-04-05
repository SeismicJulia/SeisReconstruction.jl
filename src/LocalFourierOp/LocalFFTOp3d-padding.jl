using DSP, FFTW

#=
function FFTOp(in,adj;normalize=true)
	norm = normalize ? sqrt(length(in[:])) : 1.
	if (adj)
		out = fft(in)/norm
	else
		out = bfft(in)/norm
	end

	return out
end

=#



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




function check_parameters(in,Nw,Nkx,Nky,Noverlap_t,Noverlap_x, Noverlap_y;analysis_win=nothing)
    
    
    #I should use a loop to check the parameters at once and print the right infoarmation (associated to the right variable)
    #if something fails.
    #=
     if typeof(size(in)) != Tuple{Int64, Int64}
        message="The input must be at least 2D array!"
         error(message)
     end
    =# 
    
     if Nw == 0.0
         message="Nw must be greater than 1."
         error(message)
     else
        
        Nw=Int(Nw)
        if Nw < 1
            message="Nw must be greater than 1"
            error(message)
        end
     end
  
     if Noverlap_t > Nw
        message="Noverlap in time must be less than Nw."
        error(message)
        #Noverlap = Nw ÷ 2
     else
        Noverlap_t = Int(Noverlap_t)#toma el int
        if Noverlap_t >= Nw
            message=" Noverlap must be less than Nw. The condition Noverlap < Nw is not accomplished."
            error(message)
        end
     end

     hop_t= Nw- Noverlap_t;
    
    
      if Nkx == 0.0
         message="Nkx must be greater than 1."
         error(message)
     else
        
        Nkx=Int(Nkx)
        if Nkx < 1
            message="Nk must be greater than 1"
            error(message)
        end
     end
  
     if Noverlap_x > Nkx
        message="Noverlap in space must be less than Nk."
        error(message)
        #Noverlap = Nw ÷ 2
     else
        Noverlap_x = Int(Noverlap_x)#toma el int
        if Noverlap_x >= Nkx
            message=" Noverlap must be less than Nk. The condition Noverlap < Nk is not accomplished."
            error(message)
        end
     end

     hop_x= Nkx- Noverlap_x;
    
    
    #y
    
    if Nky == 0.0
         message="Nky must be greater than 1."
         error(message)
     else
        
        Nky=Int(Nky)
        if Nky < 1
            message="Nk must be greater than 1"
            error(message)
        end
     end
  
     if Noverlap_y > Nky
        message="Noverlap in space must be less than Nky."
        error(message)
        #Noverlap = Nw ÷ 2
     else
        Noverlap_y = Int(Noverlap_y)#toma el int
        if Noverlap_y >= Nky
            message=" Noverlap must be less than Nky. The condition Noverlap < Nk is not accomplished."
            error(message)
        end
     end

     hop_y= Nky- Noverlap_y;
    
    
    
    return Int(hop_t),Int(hop_x), Int(hop_y)

end


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



function hamming3d(Nx::Int,Ny::Int,Nz::Int)
    
    x=Nx;
    y=Ny
    z=Nz;
    
    wx= hamming(x);
    wy= hamming(y);
    wz= hamming(z); 

    aux1= wx[:].* wy[:]';
    wxy=repeat(aux1 ,1,1,z ) # 64, 32 ,128
    aux2 =wz[:]'.*ones(x);
    wxz=repeat(aux2 ,1,1,y ) # 64, 32 ,128
    wxz = permutedims(wxz, [1, 3, 2]);
    
    window=wxy.*wxz;
    
    
    return window
end






function Padd3D(a,N1,N2,N3)

	n1,n2,n3 = size(a)
	b = zeros(Float64,N1,N2,N3)
	b[1:n1,1:n2,1:n3] = a

	return b
end







function stfft(in,patch_size::Tuple{Int64, Int64, Int64},Noverlap::Tuple{Int64, Int64, Int64},padd::Bool)
    
    
    Nw,Nkx,Nky= patch_size#(32,32,32)#patch_size;
    Noverlap_t, Noverlap_x, Noverlap_y = Noverlap#(16,16,16); #Number of overlap in each direction 
    hopt,hopx,hopy=check_parameters(in,Nw,Nkx,Nky,Noverlap_t,Noverlap_x, Noverlap_y; analysis_win=nothing); #Build a generic check parameters function

#windows
#wt= hamming(Nw); wx=hamming(Nkx); wy=hamming(Nky);
    #w=hamming(Nw)
    analysis_win=hamming3d(Nw,Nkx,Nky);

    Nt,Nx,Ny=size(in);
    out1=padd_boundaries(in,Nw,Nkx,Nky)
    out2,Nt2,Nx2,Ny2=padd_recover_amp(out1,Nw,Nkx,Nky,hopt,hopx,hopy)

    ω = Nw;  kx=Nkx; ky=Nky 

    τ=Int(floor((Nt2-Noverlap_t)/(Nw-Noverlap_t)));
    x=Int(floor((Nx2-Noverlap_x)/(Nkx-Noverlap_x)));
    y=Int(floor((Ny2-Noverlap_y)/(Nky-Noverlap_y))); #This is for now!


    sw=zeros(Float64,(Nw, Nkx, Nky));
    SW=zeros(ComplexF64,size(sw))
    STFT = zeros(ComplexF64,(ω,kx,ky,τ,x,y)); # preallocate the stft matrix
    
    if padd == true
            
        ωp= Nw*2;
        kxp=Nkx*2;
        kyp=Nky*2;
        STFT = zeros(ComplexF64,(ωp,kxp,kyp,τ,x,y)); # preallocate the stft matrix

        
        for k=0:y-1
            for j=0:x-1
                for i = 0:τ-1
                    sw = out2[1 + i*hopt : Nw + i*hopt, 1 + j*hopx : Nkx + j*hopx,  1 + k*hopy : Nky + k*hopy].*analysis_win;
                    sw=Padd3D(sw,ωp,kxp,kyp)
                    SW=FFTOp(sw, true, normalize=true); #normalize=false
                    STFT[ : , : , :, 1 + i , 1 + j, 1 + k]= SW[ 1 : ωp , 1 : kxp, 1 : kyp]
                end
            end
        end
        
        return STFT
        
    else
        for k=0:y-1
            for j=0:x-1
                for i = 0:τ-1
                    sw = out2[1 + i*hopt : Nw + i*hopt, 1 + j*hopx : Nkx + j*hopx,  1 + k*hopy : Nky + k*hopy].*analysis_win;
                    SW=FFTOp(sw, true, normalize=true); #normalize=false
                    STFT[ : , : , :, 1 + i , 1 + j, 1 + k]= SW[ 1 : Nw , 1 : Nkx, 1 : Nky]
                   end
                end
            end
        
        
        return STFT
    end

end






#Come back

function istfft(in,patch_size::Tuple{Int64, Int64, Int64},Noverlap::Tuple{Int64, Int64, Int64},padd::Bool; Nt=nothing,Nx=nothing, Ny=nothing)
    
    
    Nw,Nkx,Nky,=patch_size;
    Noverlap_t, Noverlap_x,Noverlap_y=Noverlap;
    
    hopt,hopx, hopy=check_parameters(in,Nw,Nkx,Nky,Noverlap_t,Noverlap_x, Noverlap_y;analysis_win=nothing)
    ω, kx,ky, τ, x, y=size(in); # size(in)
   
    # Initialize output and normalization arrays
    Ntt = Nw + (τ-1)*hopt; 
    Nxx = Nkx + (x-1)*hopx;
    Nyy = Nky + (y-1)*hopy; 
 


    synthesis_win=hamming3d(Nw,Nkx,Nky);

    aux=zeros(ComplexF64,(Ntt,Nxx,Nyy));
    norm=zeros(Float64,(Ntt,Nxx,Nyy));

    aux2=zeros(ComplexF64,size(in)) #size(in)
    xw=zeros(ComplexF64,(Nw,Nkx,Nky,τ,x,y))

     if padd == true
        
        for y=1:size(xw,6)
            for x= 1:size(xw,5)
                for τ=1:size(xw,4)
                    aux2[:,:,:,τ,x,y]=   FFTOp(in[:,:,:,τ,x,y],false,normalize=true)#/(length(in[:,:,:,τ,x,y]))
                    xw[:,:,:,τ,x,y] = aux2[1:Nw,1:Nkx,1:Nky,τ,x,y]
                end
            end
        end
        

        
        
    else
        
        for y=1:size(xw,6)
            for x= 1:size(xw,5)
                for τ=1:size(xw,4)
                    xw[:,:,:,τ,x,y]=   FFTOp(in[:,:,:,τ,x,y],false,normalize=true)
                end
            end
        end
    
    end

    
   #Re arrange eveything!

    for k=1:y
        for j=1:x
            for i=1:τ
                aux[1+(i-1)*hopt : Nw+(i-1)*hopt, 1+(j-1)*hopx : Nkx+(j-1)*hopx, 1+(k-1)*hopy : Nky+(k-1)*hopy] = 
                aux[1+(i-1)*hopt : Nw+(i-1)*hopt, 1+(j-1)*hopx : Nkx+(j-1)*hopx, 1+(k-1)*hopy : Nky+(k-1)*hopy] .+ (xw[:,:,:,i,j,k].*synthesis_win)
                norm[1+(i-1)*hopt : Nw+(i-1)*hopt, 1+(j-1)*hopx : Nkx+(j-1)*hopx, 1+(k-1)*hopy : Nky+(k-1)*hopy] = norm[1+(i-1)*hopt : Nw+(i-1)*hopt, 1+(j-1)*hopx : Nkx+(j-1)*hopx, 1+(k-1)*hopy : Nky+(k-1)*hopy] .+ (synthesis_win).^2
            end
        end
    end


    norm3=zeros(Float64,size(norm))


    for k=1:size(norm3,3)
        for j=1:size(norm3,2)
            for i=1:size(norm3,1)
                if norm[i,j,k] > 1e-10
                    norm3[i,j,k]= norm[i,j,k]
                else
                    norm3[i,j,k] = 1.0
                end
            end
        end
    end


    paddt= zeros(Float64,Nw ÷ 2)
    paddx= zeros(Float64,Nkx ÷ 2)
    paddy=zeros(Float64,Nky ÷ 2)

    aux=aux[length(paddt)+1:end-length(paddt),length(paddx)+1:end-length(paddx),length(paddy)+1:end-length(paddy)]
    norm3=norm3[length(paddt)+1:end-length(paddt),length(paddx)+1:end-length(paddx),length(paddy)+1:end-length(paddy)]


    x_out= (aux./norm3)
    x_out=x_out[1:Nt, 1:Nx, 1:Ny]
    θ=norm3[1:Nt, 1:Nx, 1:Ny]
    
    
   return x_out, θ
    
end


