g = 0

E = 70000
ρ = 2.7

F = 500
yF = 38
M = F * yF
#0-------->y
#|---------------------     "I"
#|      |     |
#V      |     |				"U"
#z      |-----|
wandstrke_U = 2
wandstrke_I = 2
breite_I = 30
breite_U = 15
hoehe_U = 20
l = 1
A_I = wandstrke_I * breite_I
A_U_l = hoehe_U * wandstrke_U
A_U_r = hoehe_U * wandstrke_U
A_U_u = breite_U * wandstrke_U - wandstrke_U * wandstrke_U * 2
A = A_I + A_U_l + A_U_r + A_U_u
y_SP_I = breite_I / 2.0
z_SP_I = wandstrke_I / 2.0
y_SP_U_l = breite_I / 2.0 - breite_U / 2.0 + wandstrke_U / 2.0
y_SP_U_r = breite_I / 2.0 + breite_U / 2.0 - wandstrke_U / 2.0
y_SP_U_u = breite_I / 2.0
z_SP_U_l = wandstrke_I + hoehe_U / 2.0
z_SP_U_r = wandstrke_I + hoehe_U / 2.0
z_SP_U_u = wandstrke_I + hoehe_U - wandstrke_U / 2.0
y_SP = (y_SP_I * A_I + y_SP_U_l * A_U_l + y_SP_U_r * A_U_r + y_SP_U_u * A_U_u) / A
z_SP = (z_SP_I * A_I + z_SP_U_l * A_U_l + z_SP_U_r * A_U_r + z_SP_U_u * A_U_u) / A
I_y_I = breite_I * wandstrke_I * wandstrke_I * wandstrke_I / 12.0
I_y_U_l = wandstrke_U * hoehe_U * hoehe_U * hoehe_U / 12.0
I_y_U_r = wandstrke_U * hoehe_U * hoehe_U * hoehe_U / 12.0
I_y_U_u = (breite_U - 2.0 * wandstrke_U) * wandstrke_U * wandstrke_U * wandstrke_U / 12.0
I_z_I = breite_I * breite_I * breite_I * wandstrke_I / 12.0
I_z_U_l = wandstrke_U * wandstrke_U * wandstrke_U * hoehe_U / 12.0
I_z_U_r = wandstrke_U * wandstrke_U * wandstrke_U * hoehe_U / 12.0
I_z_U_u = (breite_U - 2.0 * wandstrke_U) * (breite_U - 2.0 * wandstrke_U) * (breite_U - 2.0 * wandstrke_U) * wandstrke_U / 12.0
I_y = I_y_I + (z_SP_I - z_SP)^2 * A_I + I_y_U_l + (z_SP_U_l - z_SP)^2 * A_U_l + I_y_U_r + (z_SP_U_r - z_SP)^2 * A_U_r + I_y_U_u + (z_SP_U_u - z_SP)^2 * A_U_u
I_z = I_z_I + (y_SP_I - y_SP)^2 * A_I + I_z_U_l + (y_SP_U_l - y_SP)^2 * A_U_l + I_z_U_r + (y_SP_U_r - y_SP)^2 * A_U_r + I_z_U_u + (y_SP_U_u - y_SP)^2 * A_U_u
W_y = I_y / max(z_SP, wandstrke_I + hoehe_U - z_SP)

a_values = 0.00:0.01:l
a_values = [0.5]
#a_value_s = [0.25,0.33,0.5,0.66,0.75]
# maximale Biegespannung 
M_b_max = 0
a_M_b_max = 0.01

for a in a_values
	global M_b_max
	global a_M_b_max
	q = ρ * A * 1e-3 * g # *l/l
	# lagerreaktionen
	F_F_A = F * (1.0 - 3.0 * (a^2) / (l^2) + 2.0 * (a^3) / (l^3))
	F_F_B = F - F_F_A
	F_M_A = F * a * ((l - a) / l)^2
	F_M_B = F * (l - a) * (a / l)^2
	q = ρ * A * 1e-3 * g # *l/l
	q_F_A = 0.5 * q * l
	q_F_B = q_F_A
	q_M_A = -q * l^2 / 12
	q_M_B = q_M_A
	F_A = F_F_A + q_F_A
	F_B = F_F_B + q_F_B
	M_A = F_M_A + q_M_A
	M_B = F_M_B + q_M_B
	# always output abs max M 

	M_lcl =
		x -> (
			(x > a) ? (F_A * x - M_A - q * x * x / 2 - F * (x - a)) :
			((x < a) ? (F_A * x - M_A - q * x * x / 2) :
			 ([F_A * x - M_A - q * x * x / 2, F_A * x - M_A - q * x * x / 2 - F * (x - a)][argmax(abs.([F_A * x - M_A - q * x * x / 2, F_A * x - M_A - q * x * x / 2 - F * (x - a)]))])
			)
		)
	# caution: wenn mehrere x mit gleichem M_max exisiteren, wird nur eins ausgegeben
	x_canidate = 0
	if abs(M_lcl(x_canidate)) > M_b_max
		M_b_max = abs(M_lcl(x_canidate))
		a_M_b_max = a
	end
	x_canidate = l
	if abs(M_lcl(x_canidate)) > M_b_max
		M_b_max = abs(M_lcl(x_canidate))
		a_M_b_max = a
	end
	if false
		println(M_lcl(0.0))
		println(M_lcl(0.5))
		println(M_lcl(1.0))
		println(M_lcl(0.735))
		println(M_b_max)
		println(M_lcl(x_canidate))
		println("----")
	end
end

σ_b = M_b_max * 1e3 / W_y
println(a_M_b_max)

# Annamhe: Rechtecksrohr # https://www.johannes-strommer.com/formeln/flaechentraegheitsmoment-widerstandsmoment/#formel_torsion
A_t = 2.0 * (wandstrke_U * (hoehe_U + wandstrke_I)) + 2.0 * (min(wandstrke_I, wandstrke_U) * breite_U - min(wandstrke_I, wandstrke_U) * wandstrke_U * 2.0)
I_t = 2.0 * (((hoehe_U + wandstrke_I) - wandstrke_U) * (breite_U - wandstrke_U))^2.0 / (((breite_U - wandstrke_U) / wandstrke_U) + (((hoehe_U + wandstrke_I) - wandstrke_U) / wandstrke_U))
W_t = 2.0 * (breite_U - wandstrke_U) * (hoehe_U + wandstrke_I - wandstrke_U) * min((wandstrke_U), (wandstrke_I))

τ_t = M / W_t # schubspannung

σ_GEH = sqrt(σ_b^2.0 + 3.0 * τ_t^2.0) # B. Schlecht (Formel 3176)

σ_grenze = 160

SF = σ_grenze / σ_GEH
println(SF)

l=1000
durchbiegung = a->(F * a^3.0 * (l - a)^3.0) / (3.0 * l^3.0 * E * I_y)
durchbiegung_max = maximum(durchbiegung.(0:0.01:l))
durchbiegung_500 = durchbiegung.(500)
durchbiegung_667 = durchbiegung(670)
