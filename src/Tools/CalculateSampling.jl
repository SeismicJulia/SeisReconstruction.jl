"""
    CalculateSampling(in)

Calculate the sampling operator of an n-dimension input. The output has the
same size as the input. 

"""
function CalculateSampling(in)

cutoff = 1e-10
	itrace = 1

	n=size(in)
	in=reshape(in,n[1],:)
        wd = zeros(Float32,size(in))
	n2=size(in,2)
	for itrace = 1 : n2
		a = sum(in[:,itrace].*in[:,itrace])
		if (a > cutoff)
			wd[:,itrace] .= 1.
		end
	end
	wd=reshape(wd,n)
	return wd;

end
