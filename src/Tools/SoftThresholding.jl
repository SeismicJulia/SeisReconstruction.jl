
function SoftThresholding(u,η,λ)
    
    S=sign(u)*max(abs(u)- η*λ,0)
    
    return S;

end
