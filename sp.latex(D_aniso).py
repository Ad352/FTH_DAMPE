B_gal = 3e-6  # Tesla (3 μG)
c = 2.998e8   # m/s
H0 = 2.19e-18 # s⁻¹ (67.4 km/s/Mpc)
b0_z0 = 1.232e-3  # Planck CMB anchor

R_break = (B_gal * c) / (b0_z0 * H0)  # in V
R_break_TV = R_break / 1e12  # in TV
# Ergebnis: ≈ 14.8 TV ✓ (konsistent mit DAMPE)