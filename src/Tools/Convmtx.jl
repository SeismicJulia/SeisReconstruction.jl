function Convmtx(w,nr)
	w  = vec(w)
	nw = length(w)
	nd = nw + nr - 1
	cm = zeros(ComplexF64,nd, nr)
	for i=1:nr
		cm[i:i+nw-1, i] = w
	end
	return cm
end

