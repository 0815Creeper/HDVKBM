# siehe bedienhebel_festigkeit.png für eine Skizze

# |C|          _
#   ---------  D
#  |         | ¯
#  | |       |
#  | B       |
#  | |       |
#  |         |
#   ---------
#     --A--

function get_area(A, B, C, D, profile)
	if contains(profile, "H-Profil") || contains(profile, "IPE-Profil")
		return (A * D * 2 + (B - 2 * D) * C)
	elseif contains(profile, "4-Kant")
		return A * D * 2 + (B - 2 * D) * C * 2
	else
		throw(ArgumentError(string("\"", profile, "\" is not a known profile")))
	end
end

function get_widerstandsmoment(A, B, C, D, profile)
	if contains(profile, "H-Profil") || contains(profile, "IPE-Profil")
		return (A * B^3 - (A - C) * (B - 2 * D)^3) / (6 * B)
	elseif contains(profile, "4-Kant")
		return (A * B^3 - (A - 2 * C) * (B - 2 * D)^3) / (6 * B)
	else
		throw(ArgumentError(string("\"", profile, "\" is not a known profile")))
	end
end

l = 180 # siehe bedienhebel_festigkeit.png für eine Skizze
a = 32.3 # siehe bedienhebel_festigkeit.png für eine Skizze

A = get_area(20, 20, 2, 2, "4-Kant")
W = get_widerstandsmoment(20, 20, 2, 2, "4-Kant")
E = 70000
σ = 145
F = 500

function get_lager_reaktionen(F, a, l)
	# Fliegende Einzellast mit Festlager A (Hebel Drehpunkt) und Loaslager B (Anschlag unten am Schlitten)
	F_F_A = F * (l/a - 1)
	F_F_B = F - F_F_A
	F_M_A = 0
	F_M_B = 0
end

x_max_biege_M = 
M_max = 
σ_max = 
SF = 
for i in eachindex(experiments)
	global a_i = 1
	global a_max = 1
	push!(lagerreaktionen, Vector())
	push!(x_max_biege_M, zeros(length(experiments)))
	push!(M_max, zeros(length(experiments)))
	push!(σ_max, zeros(length(experiments)))
	push!(SF, zeros(length(experiments)))
	local F, ρ, A, l, W, E, σ_grenze = experiments[i]
	for a in 0:0.02:(l+0.01)
		push!(lagerreaktionen[i], get_lager_reaktionen((F, a, ρ, A, l)))
		local q = ρ * A * 1e-3 * g # *l/l
		local F_A, F_B, M_A, M_B = lagerreaktionen[i][a_i]
		# always output abs max M 
		M =
			x -> (
				(x > a) ? (F_A * x - M_A - q * x * x / 2 - F * (x - a)) :
				((x < a) ? (F_A * x - M_A - q * x * x / 2) :
				 ([F_A * x - M_A - q * x * x / 2, F_A * x - M_A - q * x * x / 2 - F * (x - a)][argmax(abs.([F_A * x - M_A - q * x * x / 2, F_A * x - M_A - q * x * x / 2 - F * (x - a)]))])
				)
			)
		# caution: wenn mehrere x mit gleichem M_max exisiteren, wird nur eins ausgegeben
		x_canidate = (F_A) / q
		if x_canidate < a && x_canidate > 0 && abs(M(x_canidate)) > abs(M(x_max_biege_M[i][a_i]))
			x_max_biege_M[i][a_i] = x_canidate
		end
		x_canidate = a
		if abs(M(x_canidate)) > abs(M(x_max_biege_M[i][a_i]))
			x_max_biege_M[i][a_i] = x_canidate
		end
		x_canidate = (F_A - F) / q
		if x_canidate > a && x_canidate < l && abs(M(x_canidate)) > abs(M(x_max_biege_M[i][a_i]))
			x_max_biege_M[i][a_i] = x_canidate
		end
		x_canidate = l
		if abs(M(x_canidate)) > abs(M(x_max_biege_M[i][a_i]))
			x_max_biege_M[i][a_i] = x_canidate
		end
		if false && a == 0.66 && F == 5*9.81 && A == 10*1*2+8*1*2
			println(A)
			println(q)
			println(F)
			println(F_A)
			println(M_A)
			println(F_B)
			println(M_B)
			println("-")
			println(M(0.0))
			println(M(0.5))
			println(M(1.0))
			println(M(0.735))
			println(x_max_biege_M[i][a_i])
			println(M(x_max_biege_M[i][a_i]))
			println("----")
		end
		if x_max_biege_M[i][a_i] < 0 || x_max_biege_M[i][a_i] > l
			throw(ArgumentError("out of bounds!"))
		end
		M_max[i][a_i] = M(x_max_biege_M[i][a_i])
		σ_max[i][a_i] = abs.(M_max[i][a_i]) * 1000 / W
		SF[i][a_i] = σ_grenze / σ_max[i][a_i]
		if SF[i][a_max] > SF[i][a_i]
			a_max = a_i
		end
		global a_i += 1
	end

	local F, ρ, A, l, W, E, σ_grenze, p1, p2, p3, p4, p5, p6 = experiments[i]
	local q = ρ * A * 1e-3 * g # *l/l
	local F_A, F_B, M_A, M_B = lagerreaktionen[i][a_max]
	x_max_biege_M[i][a_max]
	M_max[i][a_max]
	σ_max[i][a_max]
	SF[i][a_max]

	print(file, p1)
	print(file, ",")
	print(file, p2)
	print(file, ",")
	print(file, p3)
	print(file, ",")
	print(file, p4)
	print(file, ",")
	print(file, p5)
	print(file, ",")
	print(file, p6)
	print(file, ",")
	print(file, F/g)
	print(file, ",")
	print(file, collect(0:0.02:(l+0.01))[a_max])
	print(file, ",")
	print(file, σ_max[i][a_max])
	print(file, ",")
	println(file, SF[i][a_max])
end

close(file)
