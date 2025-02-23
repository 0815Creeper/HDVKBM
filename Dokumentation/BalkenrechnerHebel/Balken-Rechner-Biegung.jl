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
F = 500
wandstrke = 2
breite = 20
l = 180
A = 4 * wandstrke * breite - 4 * wandstrke * wandstrke
y_SP = breite / 2
z_SP = breite / 2
I_y = (breite^4.0 - (breite - 2 * wandstrke)^4.0) / 12.0
W_y = I_y * 2.0 / breite

x_values = 0.00:0.01:l
a = 32.3
# maximale Biegespannung 
M_b_max = 0
x_M_b_max = 0.00

# lagerreaktionen
F_A = F * (l / a - 1)
F_B = F * l / a
#M_lcl = x -> (F_A * x + F_B * (a - x) + F * (x - l))
for x in x_values
	global M_b_max
	global x_M_b_max

	M_lcl =
		x -> (
			(x < a) ? (F_A * x) :
			((x < l) ? (F_A * x + F_B * (a - x)) : 
			(F_A * x + F_B * (a - x) + F * (x - l))
			)
		)

	# caution: wenn mehrere x mit gleichem M_max exisiteren, wird nur eins ausgegeben
	if abs(M_lcl(x)) > M_b_max
		M_b_max = abs(M_lcl(x))
		x_M_b_max = x
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

σ_b = M_b_max / W_y
println(x_M_b_max)
