#=
function SoftThresholding(u,η,λ)
    
    S=sign(u)*max(abs(u)- η*λ,0)
    
    return S;

end
=#


function SoftThresholding(in,η,λ)
           
    out=sign.(in).*max.(abs.(in) .- η*λ,0)
    
    return out;

end




#=

function SoftThreshRED(u,η,λ)
    
    S=sign.(u).*max.(abs.(u) .- η*λ,0)
    
    return S;

end
=#