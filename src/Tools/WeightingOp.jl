function WeightingOp(in,adj;w=1)

	return in.*w

end

function WeightingOp(m::AbstractString,d::AbstractString,adj;w="NULL")

	if (adj==true)
		d1,h1,e1 = SeisRead(d)
		d2,h2,e2 = SeisRead(w)
		SeisWrite(m,d1[:,:].*d2[:,:],h1,e1)
	else
		d1,h1,e1 = SeisRead(m)
		d2,h2,e2 = SeisRead(w)
		SeisWrite(d,d1[:,:].*d2[:,:],h1,e1)
	end

end


