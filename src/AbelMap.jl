@inline function AbelMap_singular(x,y,w)
	return -(w/pi) * atan(2 * (pi * (x + pi * (pi / (3 * w^2)))/y) / w)
end

function GeneralAbelMap(x,y,OLS,w)
	#maybe for large x and y something else (z near zero)
	n = length(OLS)
	for i in 1:n
		q = (-(OLS[i][1] - OLS[i][2]) * 
					(OLS[i][1] - OLS[i][3]) + (x - OLS[i][1])^2/4)^0.5
		if abs(q - (x - OLS[i][1])/2) > abs(q + (x - OLS[i][1])/2)
			q = -q
		end
		x = q + (OLS[i][1] + x)/2
		y = y * (x - OLS[i][1])^2/((x - OLS[i][1])^2  - 
					(OLS[i][1] - OLS[i][2]) * (OLS[i][1] - OLS[i][3]))
	end
	return AbelMap_singular(x,y,w)
end

function CreateAbelMap(OLS,w)
	function AbelMap(x,y)
		return GeneralAbelMap(x,y,OLS,w)
	end
	return AbelMap
end
