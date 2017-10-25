function InnerProduct(in1,in2)

	return convert(Float32,real(sum(conj(in1[:]).*in2[:])))

end


