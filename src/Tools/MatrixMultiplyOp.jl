function MatrixMultiplyOp(in,adj;matrix=1)

	if (adj)
		out = matrix'*in
	else
		out = matrix*in
	end

	return out
end