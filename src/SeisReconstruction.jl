module SeisReconstruction
    using Interpolations,Requires,FFTW, LinearAlgebra,DSP
    include("Reconstruction/Reconstruction.jl")
    include("Tools/Tools.jl")
    include("Modelling/Modelling.jl")
end
