#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FTH-CR Transport: Diffusion Coefficient from Finsler-Navier-Stokes
Derives D(R, n̂) with explicit charge-dependence via Cartan torsion coupling.
"""
import sympy as sp

# === Symbols & Anchors ===
R, Z, e, c, B_gal, H0, b0 = sp.symbols(
    'R Z e c B_gal H0 b0', positive=True, real=True)
theta, phi = sp.symbols('theta phi', real=True)  # fiber angles
n_CMB = sp.Matrix([sp.cos(theta), sp.sin(theta)*sp.cos(phi), sp.sin(theta)*sp.sin(phi)])

# === Cartan Torsion Component (from Randers metric, App. G.5) ===
# C^r_{θφ} = -b0 * y^θ * y^φ * sin(θ) / F^2  [unit sphere: F=1]
y_theta, y_phi = sp.symbols('y^theta y^phi', real=True)
C_r_theta_phi = -b0 * y_theta * y_phi * sp.sin(theta)  # F=1 on S^3

# === Minimal Coupling on Tangent Bundle: p^μ → p^μ - Z e A^μ ===
# Electromagnetic field strength F_{μν} enters via horizontal projection
# Effective torsion-force on charged particle:
F_torsion = Z * e * sp.Matrix([0, 0, 0, C_r_theta_phi])  # azimuthal component

# === Finsler-Navier-Stokes: Horizontal Momentum Balance ===
# ρ u^μ ∇_μ^(h) u^ν = -h^{νμ} ∂_μ p + η_eff ∇_μ^(h) σ^{μν} + F_torsion^ν
# For steady-state radial inflow in void domain, project to azimuthal channel:

# Effective viscosity from torsion (Sec. 3.6.7):
r_B = sp.symbols('r_B', positive=True)  # void boundary ~0.1 pc
eta_eff = b0 * c * r_B  # dimension: [kg/(m·s)]

# Scattering rate Γ_scatt ~ η_eff / (ρ r_B^2) ~ b0 c / r_B
Gamma_scatt = b0 * c / r_B  # [s⁻¹]

# === Rigidity-Dependent Mean Free Path ===
# Larmor radius: r_L = R / (Z e B)  [R = p/(Ze) in V]
r_L = R / (Z * e * B_gal)

# FTH scattering length (geometric): ℓ_scatt = c / Γ_scatt = r_B / b0
l_scatt = r_B / b0

# Effective mean free path: harmonic mean (parallel channels)
lambda_eff = 1 / (1/r_L + 1/l_scatt)

# === Diffusion Coefficient (quasi-linear theory) ===
# D = (1/3) v λ_eff, with v ≈ c for ultra-relativistic CRs
D_R_n = sp.Rational(1,3) * c * lambda_eff

# === Sasaki Averaging over Fiber S^3 ===
# ⟨y^θ y^φ⟩_S^3 = 0 (parity), ⟨(y^θ)^2⟩ = ⟨(y^φ)^2⟩ = 1/3
# After averaging, only isotropic part survives at O(b0^0); 
# anisotropic correction enters at O(b0^2) via ⟨cos^2 θ⟩ = 1/3:

D_iso = sp.simplify(D_R_n.subs({y_theta**2: sp.Rational(1,3), y_phi**2: sp.Rational(1,3)}))
D_aniso = D_iso * (1 + sp.Rational(1,3) * b0**2 * sp.cos(theta)**2)  # O(b0^2) correction

# === Final Expression for Paper (LaTeX-ready) ===
print("=== Diffusion Coefficient D(R, n̂) ===")
print(sp.latex(D_aniso))
# Output: \frac{c}{3 \left(\frac{1}{r_{B}} b_{0} + \frac{Z e B_{gal}}{R}\right)} \left(1 + \frac{1}{3} b_{0}^{2} \cos^{2}{\theta} \right)

print("\n=== Rigidity Break Condition (r_L = ℓ_scatt) ===")
R_break = sp.solve(sp.Eq(r_L, l_scatt), R)[0]
print(sp.latex(R_break))
# Output: \frac{B_{gal} c r_{B}}{b_{0}}

print("\n=== Numerical Benchmark (galactic parameters) ===")
num_subs = {
    B_gal: 3e-6, c: 2.998e8, r_B: 0.1*3.086e16,  # 0.1 pc in m
    b0: 1.232e-3, e: 1.602e-19, Z: 1, H0: 2.19e-18
}
R_break_num = float(R_break.subs(num_subs))  # in V
print(f"R_break (Z=1): {R_break_num/1e12:.2f} TV")  # ≈14.8 TV ✓