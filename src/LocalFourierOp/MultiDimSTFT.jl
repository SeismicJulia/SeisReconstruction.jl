


#Updated multi dimensional windowed fourier transform

function stfft(in::Array{T, 3}, patch_size::Tuple{Int64, Int64, Int64}, Noverlap::Tuple{Int64, Int64, Int64}, padd::Bool)where T<:Number
    
    
    Nw, Nkx, Nky = patch_size
    Noverlap_t, Noverlap_x, Noverlap_y = Noverlap
    hopt, hopx, hopy = check_parameters(in, Nw, Nkx, Nky, Noverlap_t, Noverlap_x, Noverlap_y; analysis_win=nothing)

    analysis_win = hamming3d(Nw, Nkx, Nky)
    Nt, Nx, Ny = size(in)
    out1 = padd_boundaries(in, Nw, Nkx, Nky)
    out2, Nt2, Nx2, Ny2 = padd_recover_amp(out1, Nw, Nkx, Nky, hopt, hopx, hopy)

    τ = Int(floor((Nt2 - Noverlap_t) / (Nw - Noverlap_t)))
    x = Int(floor((Nx2 - Noverlap_x) / (Nkx - Noverlap_x)))
    y = Int(floor((Ny2 - Noverlap_y) / (Nky - Noverlap_y)))

    if padd
        ωp, kxp, kyp = 2 * Nw, 2 * Nkx, 2 * Nky;

        STFT = zeros(ComplexF64,(ωp, kxp, kyp, τ, x, y));
    else
        STFT = zeros(ComplexF64,(Nw, Nkx, Nky, τ, x, y));
    
    end

   @inbounds   @threads for k = 0:y-1
                       for j = 0:x-1
                             for i = 0:τ-1
                                    sw = out2[1 + i*hopt : Nw + i*hopt, 1 + j*hopx : Nkx + j*hopx, 1 + k*hopy : Nky + k*hopy] .* analysis_win
                                    if padd
                                         sw = Padd3D(sw, ωp, kxp, kyp)
                                    end
                                    SW = FFTOp(sw, true, normalize=true)
                                    if padd
                                        STFT[:, :, :, 1 + i, 1 + j, 1 + k] = SW[1:ωp, 1:kxp, 1:kyp]
                                        else
                                        STFT[:, :, :, 1 + i, 1 + j, 1 + k] = SW[1:Nw, 1:Nkx, 1:Nky]
                                    end
                            end
                        end
                    end

    return STFT
    #=
    
    Add with distributed to work in different processors.

    if padd
        ωp, kxp, kyp = 2 * Nw, 2 * Nkx, 2 * Nky
        STFT = SharedArray{ComplexF64}((ωp, kxp, kyp, τ, x, y))
    else
        STFT = SharedArray{ComplexF64}((Nw, Nkx, Nky, τ, x, y))
    end

    @distributed for k = 0:y-1
        for j = 0:x-1
            for i = 0:τ-1
                sw = out2[1 + i*hopt : Nw + i*hopt, 1 + j*hopx : Nkx + j*hopx, 1 + k*hopy : Nky + k*hopy] .* analysis_win
                if padd
                    sw = Padd3D(sw, ωp, kxp, kyp)
                end
                SW = FFTOp(sw, true, normalize=true)
                if padd
                    STFT[:, :, :, 1 + i, 1 + j, 1 + k] = SW[1:ωp, 1:kxp, 1:kyp]
                else
                    STFT[:, :, :, 1 + i, 1 + j, 1 + k] = SW[1:Nw, 1:Nkx, 1:Nky]
                end
            end
        end
    end

    return Array(STFT)

    =#
end


#Updated multi dimensional windowed fourier transform


function istfft(in::Array{ComplexF64, 6}, patch_size::Tuple{Int64, Int64, Int64}, Noverlap::Tuple{Int64, Int64, Int64}, padd::Bool; Nt=nothing, Nx=nothing, Ny=nothing)

    Nw, Nkx, Nky = patch_size
    Noverlap_t, Noverlap_x, Noverlap_y = Noverlap
    
    hopt, hopx, hopy = check_parameters(in, Nw, Nkx, Nky, Noverlap_t, Noverlap_x, Noverlap_y; analysis_win=nothing)
    ω, kx, ky, τ, x, y = size(in) # size(in)
    
    # Initialize output and normalization arrays
    Ntt = Nw + (τ-1)*hopt
    Nxx = Nkx + (x-1)*hopx
    Nyy = Nky + (y-1)*hopy
    
    synthesis_win = hamming3d(Nw, Nkx, Nky)
    
    aux = zeros(ComplexF64, (Ntt, Nxx, Nyy))
    norm = zeros(Float64, (Ntt, Nxx, Nyy))
    
    aux2 = similar(in, ComplexF64)
    xw = zeros(ComplexF64, (Nw, Nkx, Nky, τ, x, y))
    
    
    
    
    if padd
       @inbounds @threads for yi in 1:y
             for xi in 1:x
                  for ti in 1:τ
                    aux2[:, :, :, ti, xi, yi] = FFTOp(in[:, :, :, ti, xi, yi], false, normalize=true)
                    xw[:, :, :, ti, xi, yi] = aux2[1:Nw, 1:Nkx, 1:Nky, ti, xi, yi]
                end
            end
        end
    else
        
        @inbounds  @threads for yi in 1:y
                           for xi in 1:x
                               for ti in 1:τ
                                    xw[:, :, :, ti, xi, yi] = FFTOp(in[:, :, :, ti, xi, yi], false, normalize=true)
                                end
                            end
                    end
        end
    
       
    @inbounds for k=1:y
                for j=1:x
                    for i=1:τ
    
                t_index = 1 + (i-1)*hopt : Nw + (i-1)*hopt;
                x_index = 1 + (j-1)*hopx : Nkx + (j-1)*hopx;
                y_index = 1 + (k-1)*hopy : Nky + (k-1)*hopy;
                aux[t_index,x_index,y_index] = 
                aux[t_index,x_index,y_index] .+ (xw[:,:,:,i,j,k].*synthesis_win)
                norm[t_index,x_index,y_index] = norm[t_index,x_index,y_index] .+ (synthesis_win).^2
            end
        end
    end
    
    
    
    norm3 = map(x -> x > 1e-10 ? x : 1.0, norm)
    
    
    paddt = Nw ÷ 2
    paddx = Nkx ÷ 2
    paddy = Nky ÷ 2
    
    aux = aux[paddt+1:end-paddt, paddx+1:end-paddx, paddy+1:end-paddy]
    norm3 = norm3[paddt+1:end-paddt, paddx+1:end-paddx, paddy+1:end-paddy]
    
    x_out = (aux ./ norm3)
    x_out = x_out[1:Nt, 1:Nx, 1:Ny]
    θ = norm3[1:Nt, 1:Nx, 1:Ny]
    
    return x_out, θ

end

