function Periods(LS)
	d = length(LS)
	w = 1.0im * (pi / (3 * LS[d][1]) ^ (0.5))
    eta = pi*(pi / (6 * w))
    for i in 1:d
		eta = 2*eta + LS[d - i + 1][1]*w/2
	end
	return w, eta
end
