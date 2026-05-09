struct WCurve{T<:AbstractFloat,WF<:Function,SWF<:Function,RWF<:Function,AM<:Function}
	g2::Complex{T}
	g3::Complex{T}
	D::Complex{T}
	J::Complex{T}
	Roots::Tuple{Complex{T},Complex{T},Complex{T}}
	Periods::Tuple{Complex{T},Complex{T}}
	Eta::Tuple{Complex{T},Complex{T}}
	WeierstrassFunc::WF
	WeierstrassFuncSigmaSquared::SWF
	WeierstrassFuncRaw::RWF
	AbelMap::AM
	
	function WCurve(data; source = "Invariants", n = 15, e = nothing, rectangular = false)
		g2 = undef
		g3 = undef
		Roots = undef
		if source == "Invariants"
			if(length(data) != 2)
				throw(ArgumentError("Invalid length of the data argument: please, ensure that in case source==\"Invariants\", data contains (g2,g3)"))
			end
			g2 = data[1]
			g3 = data[2]
			Roots = RootsFromInvariants(g2,g3)
		elseif source == "Roots"
			if(length(data) != 3)
				throw(ArgumentError("Invalid length of the data argument: please, ensure that in case source==\"Roots\", data contains the roots (e1,e2,e3)"))
			end
			Roots = data
			s = (Roots[1] + Roots[2] + Roots[3]) + 0.0
			if(isnothing(e) && abs(s) > 10*eps(typeof(real(s))))
				throw(ArgumentError("The sum of roots does not equal to zero"))
			end
			if(!isnothing(e))
				if(abs(s) > e)
					throw(ArgumentError("The sum of roots does not equal to zero"))
				end
			end
			g2,g3 = InvariantsFromRoots(Roots)
		else
			throw(ArgumentError("Invalid source: please ensure that source is equal to either \"Invariants\", or \"Roots\";\nby default source is equal to \"Invariants\""))
		end
		###################################### only for rectangular lattices: w1 is real, w2 is purely imaginary
		#if source == "Periods"
		#	(g2,g3) = Invariants(data[1], data[2], e)
		#	Roots = roots([-g3,-g2,0,4])
		#end
		T = typeof(real(Roots[1]) + 0.0)
		if(isnothing(e))
			e = 10*eps(T)
		end
		Roots = Complex{T}.((Roots[1], Roots[2], Roots[3]))
		g2 = Complex{T}(g2)
		g3 = Complex{T}(g3)
		D,J = JInvariantAndDiscriminant(g2,g3)
		
		OLS,sOLS = ConstructLatticeSequences(Roots, e, n, T)
		
		w1,e1 = Periods(OLS)
		w2,e2 = Periods(sOLS)
		
		function WeierFuncRaw(z)
			return return WeierstrassFuncRaw(z,OLS,w1)
		end
		
		function WeierFuncSigmaSquared(z; use_periodicity = false)
			u = z
			zetashift = 0
			sigmacoeff = 1
			if (use_periodicity)
				(u,n1,n2) = shiftToFundamentalParallelogram(u,w1,w2)
				zetashift = 2*n1*e1 + 2*n2*e2
				sigmacoeff = exp(zetashift*(u + z))
			end
			wf = WeierstrassFuncSigmaSquaredAbstract(WeierFuncRaw(u))
			return (wf[1]*sigmacoeff, wf[2] + zetashift, wf[3], wf[4])
		end
		
		function WeierFunc(z; use_periodicity = false)
			u = z
			zetashift = 0
			sigmacoeff = 1
			if(use_periodicity)
				(u,n1,n2) = shiftToFundamentalParallelogram(u,w1,w2)
				zetashift = 2*n1*e1 + 2*n2*e2
				sigmacoeff = (-1)^(n1*n2 + n1 + n2)*exp(zetashift*(u + z)/2)
			end
			wf = WeierstrassFuncDuplicationAbstract(WeierFuncRaw(u/2), g2, g3)
			return (wf[1]*sigmacoeff, wf[2] + zetashift, wf[3], wf[4])
		end
		
		AbMap = CreateAbelMap(OLS, w1)
		
		u1,u2,f1,f2 = w1,w2,e1,e2
		
		if(rectangular)
			t = abs(imag(Roots[1])) + abs(imag(Roots[2])) + abs(imag(Roots[3]))
			if(abs(t) > e)
				throw(ArgumentError("The lattice is not rectangular: cannot impose rectangular normalization of periods"))
			end
			u1,u2,f1,f2 = RectangularNormalization(u1,u2,f1,f2,e)
		end
		
		return new{T, typeof(WeierFunc),typeof(WeierFuncSigmaSquared),typeof(WeierFuncRaw),typeof(AbMap)}(g2,g3,D,J,Roots,(u1,u2),(f1,f2),WeierFunc,WeierFuncSigmaSquared,WeierFuncRaw,AbMap)
	end
end

Base.show(io::IO,C::WCurve) = print(io, WCurveString(C))

function WCurveString(C::WCurve)
	return "WCurve{" * string(typeof(real(C.g2))) * "}\n" * "g2 = " * string(C.g2) * ", g3 = " * string(C.g3) * "\n" *
			"roots = " * string(C.Roots[1]) * ", " * string(C.Roots[2]) * ", " * string(C.Roots[3]) * "\n" *
			"Discriminant = " * string(C.D) * ", J-invariant = " * string(C.J) * "\n" *
			"Periods = " * string(C.Periods[1]) * ", " * string(C.Periods[2]) * "\n" *
			"Eta-periods = " * string(C.Eta[1]) * ", " * string(C.Eta[2]) * "\n" *
			"Functions: WeierstrassFunc, WeierstrassFuncSigmaSquared, WeierstrassFuncRaw, AbelMap"
end

function ConstructLatticeSequences(Roots, e, n, T)
	OLS = Array{Tuple{Complex{T},Complex{T},Complex{T}}}([])
	R = Roots
	i = 1
	while i <= n && smallest_dist(R) > e 
		R = LandenTranformOfRoots(R,choose_root_min(R))
		push!(OLS, (R[1], R[2], R[3]))
		i += 1
	end
	sOLS = Array{Tuple{Complex{T}}}([])
	R = Roots
	R = LandenTranformOfRoots(R, choose_root_mid(R))
	while (abs(R[2] - R[3]) < abs(R[1] - R[3]) && abs(R[2] - R[3]) < abs(R[1] - R[2]))
		j = choose_root_mid(R)
		push!(sOLS, (R[j],))
		R = LandenTranformOfRoots(R, j)
	end
	push!(sOLS, (R[1],))
	i = 1
	while i <= n && smallest_dist(R) > e
		R = LandenTranformOfRoots(R,choose_root_min(R))
		push!(sOLS, (R[1],))
		i += 1
	end
	return OLS, sOLS
end
