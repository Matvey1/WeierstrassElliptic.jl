@inline function WeierstrassRawSingular(z,w)
	r = pi/w
	rinv = w/pi
	ex = exp((z*r)^2/3)
	si = sin(r*z)
	co = cos(r*z)
	return rinv^2 * ex * si^2, ex*(1 - si^2/3), z*ex*si^2*2/3 + (2*rinv)*ex*si*co, ex*(1 - si^2/3)*z*r^2*2/3 - (2*r/3)*ex*si*co
end

function WeierstrassFuncRaw(z,OLS,w)
	n = length(OLS)
	S,R,Sz,Rz = WeierstrassRawSingular(z, w)
	for i in 1:n
		e1 = OLS[n-i+1][1]
		t = (OLS[n-i+1][2] - OLS[n-i+1][1])*(OLS[n-i+1][3] - OLS[n-i+1][1])
		coeff = exp(e1 * z^2)
		S,R,Sz,Rz = coeff*(R*S - e1*S^2), coeff*(R^2 - e1*R*S + t*S^2), coeff*(2*e1*z*(R*S - e1*S^2) + (Rz*S + R*Sz - 2*e1*S*Sz)), 
					coeff*(2*e1*z*(R^2 - e1*R*S + t*S^2) + (2*R*Rz - e1*(Rz*S + R*Sz) + t*2*S*Sz))
    end
	return (S,R,Sz,Rz)
end

function WeierstrassFuncSigmaSquaredAbstract(WFR) #order is (sigma^2, zeta, p, p')
	return (WFR[1], WFR[3]/(2*WFR[1]), WFR[2]/WFR[1], (WFR[4] - WFR[2]*WFR[3]/WFR[1])/WFR[1])
end

function WeierstrassFuncDuplicationAbstract(WFR, g2, g3) #we are given WFR(z/2) and need to calculate Weierstrass functions at z, order is (sigma, zeta, p, p')
	sigma = WFR[2]*WFR[3] - WFR[1]*WFR[4]
	S = sigma^2
	R = 2*WFR[2]^4 - g2*WFR[1]^2*WFR[2]^2/2 - g3*WFR[1]^3*WFR[2]/2 - (6*WFR[2]^2 - g2*WFR[1]^2/2)^2/4
	sigma = WFR[2]*WFR[3] - WFR[1]*WFR[4]
	S = sigma^2
	R = -8*WFR[2]^4 + 2*g2*WFR[1]^2*WFR[2]^2 + 2*g3*WFR[1]^3*WFR[2] + (6*WFR[2]^2 - g2*WFR[1]^2/2)^2/4
	Rz = (-32*WFR[4]*WFR[2]^3 + 4*g2*WFR[3]*WFR[1]*WFR[2]^2 + 4*g2*WFR[1]^2*WFR[4]*WFR[2] + 6*g3*WFR[1]^2*WFR[2]*WFR[3] + 2*g3*WFR[1]^3*WFR[4] +
			(6*WFR[2]^2 - g2*WFR[1]^2/2)*(12*WFR[2]*WFR[4] - g2*WFR[1]*WFR[3])/2)/2
	Sz = (4*WFR[3]*WFR[2]^3 + 12*WFR[1]*WFR[2]^2*WFR[4] - 3*g2*WFR[1]^2*WFR[3]*WFR[2] - g2*WFR[1]^3*WFR[4] - 4*g3*WFR[1]^3*WFR[3])/2
	zeta,p,pz = Sz/(2*S), R/S, (Rz - R*Sz/S)/S
	return (sigma, zeta, p, pz)
end

