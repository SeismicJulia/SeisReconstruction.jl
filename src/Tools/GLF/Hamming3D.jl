
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