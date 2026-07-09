import numpy as np
import matplotlib.pyplot as plt
import os

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../../build/iEoM/tests/bcs_result.txt')
M = np.loadtxt(filename, delimiter=" ").T

phase_a = M[0]
phase_b = M[1]
higgs_a = M[2]
higgs_b = M[3]

Delta = 3.404305e-01 # if you change the system parameters, you need to adjust Delta accordingly
E_min = (2 * Delta)**2
E_max = (2 * np.sqrt(4. + Delta**2))**2

a_theory = (E_max + E_min) / 2.
b_theory = ((E_max - E_min) / 4.)**2

omegas = np.linspace(0., 4.5, 500, dtype=complex) + 1e-4j

def terminator(z_sqr):
    z_shift = z_sqr - a_theory
    radicand = z_shift**2 - 4 * b_theory
    root = np.sqrt(np.abs(radicand), dtype=complex)
    root[z_shift.real > 2*np.sqrt(b_theory)] *= -1.
    root[np.abs(z_shift) < 2*np.sqrt(b_theory)] *= -1.j
    
    return (z_shift + root) / (2. * b_theory)

def continued_fraction(z, a, b):
    z_sqr = z**2
    r = 1. / terminator(z_sqr)
    for i in range(len(b)-1, 0, -1):
        r = z_sqr - a[i-1] - b[i] / r
    return b[0] / r

fig_term, ax_term = plt.subplots()
T = terminator(omegas**2)
ax_term.plot(omegas.real, T.real, label=r"$\Re T(\omega)$")
ax_term.plot(omegas.real, T.imag, label=r"$\Im T(\omega)$")
ax_term.legend()
ax_term.set_xlabel(r"$\omega$")
ax_term.set_ylabel(r"$T(\omega)$")

###############################

fig_coeff, ax_coeff = plt.subplots(nrows=2, sharex=True)

ax_coeff[0].plot(phase_a / a_theory, "o-", label="Phase")
ax_coeff[0].plot(higgs_a / a_theory, "^-", label="Higgs")
ax_coeff[0].axhline(1, ls="--", c="k", label=r"$\infty$")

ax_coeff[1].plot(phase_b / b_theory, "o-", label="Phase")
ax_coeff[1].plot(higgs_b / b_theory, "^-", label="Higgs")
ax_coeff[1].axhline(1, ls="--", c="k", label=r"$\infty$")

ax_coeff[0].set_ylim(0.7, 1.3)
ax_coeff[1].set_ylim(0.7, 1.3)
ax_coeff[1].set_xlabel("Iteration $i$")
ax_coeff[0].set_ylabel(r"$a_i / a_\infty$")
ax_coeff[1].set_ylabel(r"$b_i^2 / b_\infty^2$")

ax_coeff[0].legend(loc="upper right")

###############################

fig, ax = plt.subplots()
g_phase = continued_fraction(omegas, phase_a, phase_b)
g_higgs = continued_fraction(omegas, higgs_a, higgs_b)

ax.plot(omegas.real, -g_phase.imag, label="Phase")
ax.plot(omegas.real, -g_higgs.imag, label="Higgs")

ax.text(0.15, 4., "Goldstone mode", color="C0", rotation=90, va="top", ha="left")
ax.text(0.9, 4.2, "Higgs Mode", color="C1", rotation=90, va="top", ha="left")

ax.set_ylim(0, 5)
ax.set_xlabel(r"$\omega$")
ax.set_ylabel(r"$A(\omega)$")
ax.legend()

plt.show()