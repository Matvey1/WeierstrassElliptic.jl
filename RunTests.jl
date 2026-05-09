module WeierstrassEllipticTests

using WeierstrassElliptic
using Test
using Random
using EllipticFunctions:wp, wsigma
import SpecialFunctions:gamma

@testset "Input correctness" begin
	@test_throws ArgumentError WCurve([1,2,3])
	@test_throws ArgumentError WCurve([1])
	@test_throws ArgumentError WCurve([1,2,3], source = "Roots")
	@test_throws ArgumentError WCurve([1,2], source = "Roots")
	@test_throws ArgumentError WCurve([1,2,3,4], source = "Roots")
	@test_throws ArgumentError WCurve([1,-1-1*im,1*im], source = "Roots", rectangular = true)
	@test_throws ArgumentError WCurve([1,2,3], source = "something")
end

@testset "Basic tests of periods" begin
	#A computable example
	C = WCurve([-1,0,1], source = "Roots", rectangular = true)
	@test abs(C.Periods[1] - gamma(1/4)^2/(2*sqrt(2*pi))) < 1e-14
	@test abs(C.Periods[2]/im - gamma(1/4)^2/(2*sqrt(2*pi))) < 1e-14
	#Legendre's identity
	r1 = rand() + im*rand()
	r2 = rand() + im*rand()
	r3 = -r1 - r2
	C = WCurve([r1,r2,r3], source = "Roots")
	@test abs(C.Periods[2] * C.Eta[1] - C.Periods[1]*C.Eta[2] - sign(imag(C.Periods[2]/C.Periods[1]))*im*pi) < 1e-12
	setprecision(256)
	r1 = BigFloat(rand()) + im*rand()
	r2 = rand() + im*rand()
	r3 = -r1 - r2
	C = WCurve([r1,r2,r3], source = "Roots")
	@test abs(C.Periods[2] * C.Eta[1] - C.Periods[1]*C.Eta[2] - sign(imag(C.Periods[2]/C.Periods[1]))*im*pi) < 1e-50
end

@testset "Basic identities for Weierstrass functions" begin
	r1 = rand() + 0.1
	r2 = im*(rand() + 0.1)
	r3 = -r1 - r2
	C = WCurve([r1,r2,r3], source = "Roots")
	z = -(0.25 + rand()/2)*C.Periods[1]/2 + -(0.25 + rand()/2)*C.Periods[2]/2
	(sigma, zeta,p,pz) = C.WeierstrassFunc(z)
	#coincidence with EllipticFunctions.jl realization
	@test abs(sigma - wsigma(z, g = (C.g2, C.g3))) < 1e-12
	@test abs(p - wp(z, g = (C.g2, C.g3))) < 1e-12
	#Differential equation of Weierstrass P
	@test abs(pz^2 - 4 * p^3 + C.g2 * p + C.g3) < 1e-12
	#Periodicity tests
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[1] + sigma*exp(C.Eta[1]*(2*z + C.Periods[1]))) < 1e-12
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[2] - zeta - 2*C.Eta[1])/abs(zeta) < 1e-12
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[3] - p)/abs(p) < 1e-12
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[4] - pz)/abs(pz) < 1e-12
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[1] + sigma*exp(C.Eta[2]*(2*z + C.Periods[2]))) < 1e-12
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[2] - zeta - 2*C.Eta[2])/abs(zeta) < 1e-12
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[3] - p)/abs(p) < 1e-12
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[4] - pz)/abs(pz) < 1e-12
	#use_periodicity test
	w = z + C.Periods[1] - C.Periods[2]
	WF = C.WeierstrassFunc(w)
	WFalt = C.WeierstrassFunc(w, use_periodicity = true)
	@test sum([abs(WF[j] - WFalt[j])/abs(WF[j]) for j in 1:4])/4 < 1e-12
	
	#tests for BigFloat compatibility
	setprecision(256)
	r1 = rand() + BigFloat(0.1)
	r2 = im*(rand() + 0.1)
	r3 = -r1 - r2
	C = WCurve([r1,r2,r3], source = "Roots")
	z = -(0.25 + rand()/2)*C.Periods[1]/2 + -(0.25 + rand()/2)*C.Periods[2]/2
	(sigma, zeta,p,pz) = C.WeierstrassFunc(z)
	@test abs(pz^2 - 4 * p^3 + C.g2 * p + C.g3) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[1] + sigma*exp(C.Eta[1]*(2*z + C.Periods[1]))) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[2] - zeta - 2*C.Eta[1])/abs(zeta) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[3] - p)/abs(p) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[1])[4] - pz)/abs(pz) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[1] + sigma*exp(C.Eta[2]*(2*z + C.Periods[2]))) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[2] - zeta - 2*C.Eta[2])/abs(zeta) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[3] - p)/abs(p) < 1e-50
	@test abs(C.WeierstrassFunc(z + C.Periods[2])[4] - pz)/abs(pz) < 1e-50
	w = z + C.Periods[1] - C.Periods[2]
	WF = C.WeierstrassFunc(w)
	WFalt = C.WeierstrassFunc(w, use_periodicity = true)
	@test sum([abs(WF[j] - WFalt[j])/abs(WF[j]) for j in 1:4])/4 < 1e-50
end

@testset "Abel map test" begin
	#test that Ab is the inverse of wp
	for i in 1:5
		r1 = rand() + 0.1
		r2 = im*(rand() + 0.1)
		r3 = -r1 - r2
		C = WCurve([r1,r2,r3], source = "Roots")
		z = (0.25 + rand()/2)*C.Periods[1]/2 + (0.25 + rand()/2)*C.Periods[2]/2
		WF = C.WeierstrassFunc(z)
		x,y = WF[3],WF[4]
		w = C.AbelMap(x,y)-z
		t1 = imag(w*conj(C.Periods[2]))/imag(C.Periods[1] * conj(C.Periods[2]))
		t2 = -imag(w*conj(C.Periods[1]))/imag(C.Periods[1] * conj(C.Periods[2]))
		d1 = t1 - floor(t1)
		d2 = t2 - floor(t2)
		if 2 * d1 > 1
			d1 = d1 - 1
		end
		if 2 * d2 > 1
			d2 = d2 - 1
		end
		u = d1*C.Periods[1] + d2*C.Periods[2]
		@test abs(u) < 1e-12
	end
	#test for BigFloat compatibility
	setprecision(256)
	for i in 1:5
		r1 = rand() + BigFloat(0.1)
		r2 = im*(rand() + 0.1)
		r3 = -r1 - r2
		C = WCurve([r1,r2,r3], source = "Roots")
		z = (0.25 + rand()/2)*C.Periods[1]/2 + (0.25 + rand()/2)*C.Periods[2]/2
		WF = C.WeierstrassFunc(z)
		x,y = WF[3],WF[4]
		w = C.AbelMap(x,y)-z
		t1 = imag(w*conj(C.Periods[2]))/imag(C.Periods[1] * conj(C.Periods[2]))
		t2 = -imag(w*conj(C.Periods[1]))/imag(C.Periods[1] * conj(C.Periods[2]))
		d1 = t1 - floor(t1)
		d2 = t2 - floor(t2)
		if 2 * d1 > 1
			d1 = d1 - 1
		end
		if 2 * d2 > 1
			d2 = d2 - 1
		end
		u = d1*C.Periods[1] + d2*C.Periods[2]
		@test abs(u) < 1e-50
	end
end

end #module