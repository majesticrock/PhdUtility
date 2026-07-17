import numpy as np
import matplotlib.pyplot as plt
import os

import matplotlib as mpl
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.pad_inches'] = 0.01
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

FIG_X = 5
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

#fig_term, ax_term = plt.subplots()
#T = terminator(omegas**2)
#ax_term.plot(omegas.real, T.real, label=r"$\Re T(\omega)$")
#ax_term.plot(omegas.real, T.imag, label=r"$\Im T(\omega)$")
#ax_term.legend()
#ax_term.set_xlabel(r"$\omega$")
#ax_term.set_ylabel(r"$T(\omega)$")

###############################

fig_coeff, ax_coeff = plt.subplots(nrows=2, sharex=True, figsize=(FIG_X, 0.75*FIG_X), layout="constrained")

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

fig, ax = plt.subplots(figsize=(FIG_X, 0.75*FIG_X), layout="constrained")
g_phase = continued_fraction(omegas, phase_a, phase_b)
g_higgs = continued_fraction(omegas, higgs_a, higgs_b)

ax.plot(omegas.real, -g_phase.imag / np.pi, label="Phase")
ax.plot(omegas.real, -g_higgs.imag / np.pi, label="Higgs")

ax.text(0.15, 1.5, "Goldstone mode", color="C0", rotation=90, va="top", ha="left")
ax.text(0.9, 1.7, "Higgs Mode", color="C1", rotation=90, va="top", ha="left")

ax.set_ylim(0, 2)
ax.set_xlabel(r"$\omega$")
ax.set_ylabel(r"$A(\omega)$")
ax.legend()

###############################

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../../build/iEoM/tests/bcs_goldstone_result.txt')
with open(filename, "r") as f:
    ED_DATA = [list(map(float, line.split())) for line in f]
EVS         = np.array(ED_DATA[0])
WEIGHTS     = np.array(ED_DATA[1])
AMPLITUDES  = np.array(ED_DATA[2])

fig_ed, ax_ed = plt.subplots(ncols=2, figsize=(2*FIG_X, 0.75*FIG_X), layout="constrained")
ax_ed[0].set_xlabel(r"$\omega$")
ax_ed[0].set_ylabel(r"$A_\mathrm{Phase} (\omega)$")

ax_ed[1].set_xlabel(r"$k / \pi$")
ax_ed[1].set_ylabel(r"$\psi_k / \max |\psi_k|$")

omegas = np.linspace(0., 4.5, 5000, dtype=complex) + 3e-2j
ks = np.linspace(-1, 1, len(AMPLITUDES))
spectral_phase = np.array([ np.sum(WEIGHTS / (z**2 - EVS**2)) for z in omegas])

ax_ed[0].plot(omegas.real, -spectral_phase.imag / np.pi)
ax_ed[1].plot(ks, AMPLITUDES / np.max(np.abs(AMPLITUDES)))
ax_ed[0].set_ylim(0, 2)
ax_ed[1].set_xticks([-1, -0.5, 0, 0.5, 1])
ax_ed[1].axvline(0.5, c="k", ls=":")
ax_ed[1].axvline(-0.5, c="k", ls=":")

fig.savefig(os.path.join(dirname, "spectral.pdf"))
fig_coeff.savefig(os.path.join(dirname, "lanczos.pdf"))
fig_ed.savefig(os.path.join(dirname, "amplitudes.pdf"))

plt.show()