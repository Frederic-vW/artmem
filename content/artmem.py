"""
Artificial membrane practical: estimate the relative permeability P_K/P_Cl
from measured potential differences across a KCl-selective membrane.

Setup: two chambers separated by an ion-selective membrane.
  - ICF chamber: KCl concentration fixed at 150 mmol/L.
  - ECF chamber: KCl concentration varied across a fixed set of values.
Students record the potential difference (mV) between the chambers for
each ECF concentration they test, then fit the Nernst-Planck relationship
to estimate P_K/P_Cl.

Reference for tabulated KCl activities: Hamer & Wu (1972), "Osmotic
Coefficients and Mean Activity Coefficients of Uni-univalent Electrolytes
in Water at 25 degC", J. Phys. Chem. Ref. Data 1(4):1047 (as used by the
original ArtMem.exe tool). Values are expressed directly in mmol/L, which
is an excellent approximation to mmol/kg for these dilute solutions.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- physical constants -----------------------------------------------
R = 8.3145  # J / (mol K)
F = 96485.0  # C / mol
T0 = 273.15  # 0 degC in K
T = 20 + T0  # fixed room temperature for this practical (20 degC)
cRTF = R * T / F / np.log(10)  # V, RT/F expressed for use with log10

# --- practical setup ----------------------------------------------------
ICF_CONCENTRATION = 150  # mmol/L, fixed for this practical
ECF_CONCENTRATIONS = [1.5, 5, 10, 50, 100, 150]  # mmol/L, the choices students test

# --- tabulated KCl activities (mmol/L), see module docstring -----------
ACTIVITY_TABLE = {
    1: 0.9670,
    1.5: 1.4408,
    5: 4.6500,
    10: 9.0500,
    50: 41.0500,
    100: 77.3000,
    150: 111.9750,
}


def davies_activity(c_mM, A=0.509):
    """
    Estimate KCl activity (mmol/L) at 25 degC from the Davies equation, a
    Debye-Hueckel extension valid up to ionic strength ~0.5 mol/L.

    Used as a fallback for concentrations not in ACTIVITY_TABLE. Matches
    the tabulated values to within ~2% over 1.5-150 mmol/L (worst at the
    high end, as expected for the Davies approximation).
    """
    I = c_mM * 1e-3  # mol/L; ionic strength of a 1:1 electrolyte = its molarity
    log10_gamma = -A * (np.sqrt(I) / (1 + np.sqrt(I)) - 0.3 * I)
    gamma = 10 ** log10_gamma
    return c_mM * gamma


def get_activity(c_mM):
    """KCl activity (mmol/L) for concentration c_mM: tabulated value if
    available, otherwise a Davies-equation estimate."""
    if c_mM in ACTIVITY_TABLE:
        return ACTIVITY_TABLE[c_mM]
    return davies_activity(c_mM)


def nernst_potential_mV(c_out, c_in, z=1):
    """Nernst equilibrium potential (mV) for an ion of valence z."""
    return 1e3 * cRTF / z * np.log10(np.asarray(c_out) / c_in)


def fit_permeability_ratio(KCl_ecf, V_exp, KCl_icf=ICF_CONCENTRATION, r_true=None, ax=None):
    """
    Estimate P_K/P_Cl from measured potentials.

    KCl_ecf : array-like of ECF KCl concentrations (mmol/L)
    V_exp   : array-like of measured potentials (mV), same length as KCl_ecf
    KCl_icf : ICF KCl concentration (mmol/L), fixed for this practical
    r_true  : optional known r, drawn on the plot title for comparison
              (only meaningful for simulated/demo data)
    ax      : optional matplotlib Axes to draw into (a new figure is
              created if None)

    Returns a dict with r_est, slope_fit_mV_per_decade, V_fit, q_a, ax.
    """
    KCl_ecf = np.asarray(KCl_ecf, dtype=float)
    V_exp = np.asarray(V_exp, dtype=float)
    if len(KCl_ecf) < 2:
        raise ValueError("Need at least two data points to fit a slope.")

    a_ecf = np.array([get_activity(c) for c in KCl_ecf])
    a_icf = get_activity(KCl_icf)
    q_a = a_ecf / a_icf

    # V = slope * log10(q_a) + offset;  theory: slope = (r-1)/(r+1) * cRTF
    slope_permV, offset_mV = np.polyfit(np.log10(q_a), V_exp, deg=1)
    slope_V = 1e-3 * slope_permV
    r_est = (1 + slope_V / cRTF) / (1 - slope_V / cRTF)
    V_fit = slope_permV * np.log10(q_a) + offset_mV

    # pick the Nernst limit curve (pure-K or pure-Cl electrode) matching
    # the sign of the fitted permeability ratio, for reference on the plot
    z = 1 if r_est > 1 else -1
    V_nernst = nernst_potential_mV(a_ecf, a_icf, z=z)

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    order = np.argsort(q_a)
    label_nernst = "Nernst pot. " + (r"$K^{+}$" if z == 1 else r"$Cl^{-}$")
    ax.semilogx(q_a[order], V_nernst[order], ':ok', mfc='none', label=label_nernst)
    ax.semilogx(q_a, V_exp, 'sb', mfc='b', label=r"measured $V_{exp}$")
    ax.semilogx(q_a[order], V_fit[order], '-dr', mfc='none', ms=8, label='fit')
    ax.set_xlabel(r"$a[KCl]_{ecf} \, / \, a[KCl]_{icf}$")
    ax.set_ylabel("V (mV)")
    ax.legend(facecolor=(0.8, 0.8, 0.8, 0.5), shadow=True, framealpha=1)
    ax.grid(visible=True, which='both', alpha=0.3)
    r_str = f"{r_true:.2f}" if r_true else "unknown"
    ax.set_title(r"Estimated $P_K/P_{Cl}$" + f" = {r_est:.2f} (true: {r_str})")

    return {
        "r_est": r_est,
        "slope_fit_mV_per_decade": slope_permV,
        "V_fit": V_fit,
        "q_a": q_a,
        "ax": ax,
    }


def plot_KCl_activity(concentrations=ECF_CONCENTRATIONS, ax=None):
    """Plot tabulated/estimated KCl activity vs. concentration, for comparison."""
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 3))
    concentrations = sorted(concentrations)
    ax.plot(concentrations, [get_activity(c) for c in concentrations],
             '-ok', mfc='none', label='KCl activity')
    ax.plot(concentrations, concentrations, ':k', label='KCl conc.')
    ax.set_xlabel("KCl conc. (mmol/L)")
    ax.set_ylabel("KCl act. (mmol/L)")
    ax.legend(loc='upper left')
    return ax
