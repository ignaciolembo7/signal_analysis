#NMRSI - Ignacio Lembo Ferrari - 02/09/2024

import numpy as np

_GYRO = 267.52218744  # gyromagnetic ratio of 1H, units: ms^-1 mT^-1; [D0] = m2/ms

# ---------------------------------------------------------------------------
# Anomalous diffusion signal model  
# ---------------------------------------------------------------------------

def M_anom_ogse(TE, G, N, M0, beta, C):
    return M0*np.exp((-1.0/12.0) * (_GYRO**2) * (G**2) * C * TE**(beta+1) / ((N)**beta))


# ---------------------------------------------------------------------------
# PGSE signal models
# ---------------------------------------------------------------------------

def M_pgse_free(TE, G, M0, alpha, D0):
    return M0*np.exp( (-1.0/12.0) * _GYRO**2 * G**2 * (alpha*D0) * TE**3)

def M_pgse_rest(TE, G, tc, M0, D0):
    return M0*np.exp( (-1.0) * _GYRO**2 * G**2 * (D0*(tc**2)*TE - D0*(tc**3)*(1 - np.exp(-TE/(2*tc)))**2))

# ---------------------------------------------------------------------------
# Internal helper: restricted-diffusion log-attenuation terms
# ---------------------------------------------------------------------------

def _rest_log_attenuation(N, x, y, tc, bSE):
    """Return the three additive log-attenuation terms for the Callaghan
    restricted-diffusion model shared by both NOGSE and OGSE signal functions.

    Parameters are already cast to numpy arrays by the caller.
    The signal is:  M0 * exp(-phi_SE) * exp(-phi_N) * exp(phi_cross)

    Physical meaning of the three terms:
      phi_SE    -- SE lobe (last lobe, duration y)
      phi_N     -- N-1 NOGSE/OGSE lobes
      phi_cross -- cross-correlation between the two groups of lobes
    """
    # N=1: no oscillating lobes — only the SE lobe contributes (phi_N = phi_cross = 0)
    if int(N) == 1:
        phi_SE = bSE**2 * tc**2 * (
            4 * np.exp(-y / tc / 2)
            - np.exp(-y / tc)
            - 3
            + y / tc
        )
        return phi_SE, np.zeros_like(phi_SE), np.zeros_like(phi_SE)

    # --- sub-expressions that repeat inside the NOGSE/OGSE lobe sum ---
    e_train     = np.exp(-(N - 1) * x / tc)          # decay over all N-1 lobes
    e_lobe      = e_train ** (1 / (N - 1))            # single-lobe decay  = exp(-x/tc)
    e_lobe_half = e_train ** (1 / (N - 1) / 2)        # half-lobe decay    = exp(-x/(2*tc))
    neg_lobe_N  = (-e_lobe) ** (N - 1)                # sign-alternating train factor
    D           = e_lobe + 1                           # denominator in interference sums

    # --- phase from the N-1 NOGSE/OGSE lobes ---
    phi_N = bSE**2 * tc**2 * (
        (N - 1) * x / tc
        + (-1) ** (N - 1) * e_train
        + 1 - 2 * N
        - 4 * e_lobe_half * neg_lobe_N     / D
        + 4 * e_lobe_half                  / D
        + 4 * neg_lobe_N  * e_lobe         / D**2
        + 4 * e_lobe * ((N - 1) * e_lobe + N - 2) / D**2
    )

    # --- phase from the final SE lobe (duration y) ---
    phi_SE = bSE**2 * tc**2 * (
        4 * np.exp(-y / tc / 2)
        - np.exp(-y / tc)
        - 3
        + y / tc
    )

    # --- cross-correlation between NOGSE lobes and SE lobe ---
    inner_cross = (
        (
            np.exp((-y + 2 * x) / tc / 2)
            + np.exp((x - 2 * y) / tc / 2)
            - np.exp((x - y) / tc) / 2
            - np.exp(-y / tc) / 2
            + np.exp(x / tc / 2)
            + np.exp(-y / tc / 2)
            - np.exp(x / tc) / 2
            - 0.1e1 / 0.2e1
        ) * (-1) ** (2 * N)
        + 2 * (-1) ** (1 + N) * np.exp(-(2 * N * x - 3 * x + y) / tc / 2)
        + (
            np.exp(((3 - 2 * N) * x - 2 * y) / tc / 2)
            - np.exp((-N * x + 2 * x - y) / tc) / 2
            + np.exp(-(2 * N * x - 4 * x + y) / tc / 2)
            + np.exp(-(2 * N * x - 2 * x + y) / tc / 2)
            - np.exp((-N * x + x - y) / tc) / 2
            + np.exp(-x * (-3 + 2 * N) / tc / 2)
            - np.exp(-x * (N - 2) / tc) / 2
            - e_train / 2
        ) * (-1) ** N
        + 2 * (-1) ** (1 + 2 * N) * np.exp((x - y) / tc / 2)
    )
    phi_cross = 2 * tc**2 * inner_cross * bSE**2 / (np.exp(x / tc) + 1)

    return phi_SE, phi_N, phi_cross


# ---------------------------------------------------------------------------
# NOGSE signal models
# ---------------------------------------------------------------------------

def M_nogse_free(TE, G, N, x, M0, D0):
    x  = np.array(x)
    TE = np.array(TE)
    N  = np.array(N)
    G  = np.array(G)

    y = TE - (N - 1) * x

    return M0 * np.exp( (-1.0 / 12) * _GYRO**2 * G**2 * D0 * ((N - 1) * x**3 + y**3))


def M_nogse_free_offset(TE, G, N, x, M0, D0, C):
    return M_nogse_free(TE, G, N, x, M0, D0) + C


def M_nogse_rest(TE, G, N, x, tc, M0, D0):
    x  = np.array(x)
    TE = np.array(TE)
    N  = np.array(N)
    G  = np.array(G)

    y   = TE - (N - 1) * x
    bSE = _GYRO * G * np.sqrt(D0 * tc)

    phi_SE, phi_N, phi_cross = _rest_log_attenuation(N, x, y, tc, bSE)
    return M0 * np.exp(-phi_SE) * np.exp(-phi_N) * np.exp(phi_cross)


def M_nogse_rest_offset(TE, G, N, x, tc, M0, D0, C):
    return C + M_nogse_rest(TE, G, N, x, tc, M0, D0)


def M_nogse_mixed(TE, G, N, x, tc, alpha, M0, D0):  # alpha is 1/alpha
    return M0 * M_nogse_free(TE, G, N, x, 1, alpha * D0) * M_nogse_rest(TE, G, N, x, tc, 1, (1 - alpha) * D0)


def M_nogse_mixto_offset(TE, G, N, x, tc, alpha, M0, D0, C):  # alpha is 1/alpha
    return M0 * M_nogse_mixed(TE, G, N, x, tc, alpha, 1, D0) + C


# ---------------------------------------------------------------------------
# NOGSE contrast models
# ---------------------------------------------------------------------------

def NOGSE_contrast_vs_g_free(TE, G, N, M0, D0):
    return M_nogse_free(TE, G, N, TE / N, M0, D0) - M_nogse_free(TE, G, N, 0, M0, D0)


def NOGSE_contrast_vs_g_free_grad_offset(TE, G, N, M0, D0, g0_mTm):
    G_eff = np.array(G) + float(g0_mTm)
    return NOGSE_contrast_vs_g_free(TE, G_eff, N, M0, D0)


def NOGSE_contrast_vs_g_rest(TE, G, N, tc, M0, D0):
    return M_nogse_rest(TE, G, N, TE / N, tc, M0, D0) - M_nogse_rest(TE, G, N, 0, tc, M0, D0)


def NOGSE_contrast_vs_g_tort(TE, G, N, alpha, M0, D0):  # alpha is 1/alpha
    return M_nogse_free(TE, G, N, TE / N, M0, alpha * D0) - M_nogse_free(TE, G, N, 0, M0, alpha * D0)


def NOGSE_contrast_vs_g_mixed(TE, G, N, tc, alpha, M0, D0):
    return M_nogse_mixed(TE, G, N, TE / N, tc, alpha, M0, D0) - M_nogse_mixed(TE, G, N, 0, tc, alpha, M0, D0)


def NOGSE_contrast_vs_ad(Lc, Ld, n, alpha, D0):
    # Lc and Ld are inverted, and alpha is 1/alpha.
    # Mathematica-generated closed-form expression for the NOGSE contrast
    # parametrized by compartment length scales instead of tc.

    # Argument of the reference exponential (NOGSE x=0 limit):
    exponent_ref = D0**3 * (
        (-0.08333333333333333 * alpha * Lc**6) / D0**3
        - (
            (1 - alpha) * Ld**4
            * (
                Lc**2 / D0
                + ((-3 - np.e**(-Lc**2 / Ld**2) + 4 / np.e**(Lc**2 / (2.0 * Ld**2))) * Ld**2) / D0
            )
        ) / D0**2
    )

    # Argument of the NOGSE exponential (NOGSE x=TE/N limit):
    exponent_nogse = D0**3 * (
        (
            2 * (-1)**n * (1 - alpha) * (
                -3.0 * (-1)**n
                - 1 / (2.0 * np.e**(Lc**2 / Ld**2))
                - (0.5 * (-1)**n) / np.e**(Lc**2 / (Ld**2 * n))
                + (2.0 * (-1)**n) / np.e**(Lc**2 / (2.0 * Ld**2 * n))
                + 2.0 * (-1)**n * np.e**(Lc**2 / (2.0 * Ld**2 * n))
                - 0.5 * (-1)**n * np.e**(Lc**2 / (Ld**2 * n))
                + 2 * np.e**((Lc**2 * (3 - 2 * n)) / (2.0 * Ld**2 * n))
                - 1 / (2.0 * np.e**((Lc**2 * (-2 + n)) / (Ld**2 * n)))
                + 2 * np.e**((D0 * (Lc**2 / D0 - (2 * Lc**2 * n) / D0)) / (2.0 * Ld**2 * n))
                - 3 * np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n))
            ) * Ld**6
        ) / (D0**3 * (1 + np.e**(Lc**2 / (Ld**2 * n))))
        - (0.08333333333333333 * alpha * Lc**6) / (D0**3 * n**2)
        - (
            (-1 + alpha) * Ld**4
            * (
                -((np.e**(Lc**2 / (Ld**2 * n)) * Lc**2) / D0)
                + ((1 - 4 * np.e**(Lc**2 / (2.0 * Ld**2 * n)) + 3 * np.e**(Lc**2 / (Ld**2 * n))) * Ld**2 * n) / D0
            )
        ) / (D0**2 * np.e**(Lc**2 / (Ld**2 * n)) * n)
        - (
            (1 - alpha) * Ld**6
            * (
                1
                + (-1)**(1 + n) * np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n))
                - (
                    4 * (-(np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n)))**n
                ) / (1 + (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n)))**2
                + (
                    4 * (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (2.0 * (-1 + n)))
                ) / (1 + (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n)))
                + (
                    4
                    * (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (2 - 2 * n))
                    * (-(np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n)))**n
                ) / (1 + (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n)))
                + (Lc**2 * (-1 + n)) / (Ld**2 * n)
                - 2 * n
                + (
                    4
                    * (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n))
                    * (
                        -2
                        + (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n)) * (-1 + n)
                        + n
                    )
                ) / (1 + (np.e**((D0 * (Lc**2 / D0 - (Lc**2 * n) / D0)) / (Ld**2 * n)))**(1 / (-1 + n)))**2
            )
        ) / D0**3
    )

    return -np.exp(exponent_ref) + np.exp(exponent_nogse)


# ---------------------------------------------------------------------------
# OGSE signal models
# ---------------------------------------------------------------------------

def M_ogse_free(TE, G, N, x, M0, D0):
    x  = np.array(x)
    TE = np.array(TE)
    N  = np.array(N)
    G  = np.array(G)

    y = TE - (N - 1) * x

    return M0 * np.exp(-1.0 / 12 * _GYRO**2 * G**2 * D0 * ((N - 1) * x**3 + y**3))

def M_ogse_tort(TE, G, N, x, alpha, M0, D0):
    return M_ogse_free(TE, G, N, x, M0, alpha*D0)

def M_ogse_rest(TE, G, N, x, tc, M0, D0):
    x  = np.array(x)
    TE = np.array(TE)
    N  = np.array(N)
    G  = np.array(G)

    y   = TE - (N - 1) * x
    bSE = _GYRO * G * np.sqrt(D0 * tc)

    phi_SE, phi_N, phi_cross = _rest_log_attenuation(N, x, y, tc, bSE)
    return M0 * np.exp(-phi_SE) * np.exp(-phi_N) * np.exp(phi_cross)


def M_ogse_rest_offset(TE, G, N, x, tc, M0, D0, C):
    clean_signal = M_ogse_rest(TE, G, N, x, tc, M0, D0)
    return np.sqrt(clean_signal**2 + C**2) / np.sqrt(1 + C**2)


def M_ogse_mixed(TE, G, N, x, tc, alpha, M0, D0):  # alpha is 1/alpha
    return M0 * M_nogse_free(TE, G, N, x, 1, alpha * D0) * M_nogse_rest(TE, G, N, x, tc, 1, (1 - alpha) * D0)


def M_ogse_mixed_offset(TE, G, N, x, tc, alpha, M0, D0, C, RN):
    clean_signal = M_ogse_mixed(TE, G, N, x, tc, alpha, M0, D0) + C
    return np.sqrt(clean_signal**2 + RN**2)


# ---------------------------------------------------------------------------
# OGSE contrast models
# ---------------------------------------------------------------------------

def OGSE_contrast_vs_g_free(TE, G1, G2, N1, N2, M0, D0):
    return M_ogse_free(TE, G1, N1, TE / N1, M0, D0) - M_ogse_free(TE, G2, N2, TE / N2, M0, D0)


def OGSE_contrast_vs_g_rest(TE, G1, G2, N1, N2, tc, M0, D0):
    return M_ogse_rest(TE, G1, N1, TE / N1, tc, M0, D0) - M_ogse_rest(TE, G2, N2, TE / N2, tc, M0, D0)


def OGSE_contrast_vs_g_rest_offset(TE, G1, G2, N1, N2, tc, M0, D0, C):
    return M_ogse_rest_offset(TE, G1, N1, TE / N1, tc, M0, D0, C) - M_ogse_rest_offset(TE, G2, N2, TE / N2, tc, M0, D0, C)


def OGSE_contrast_vs_g_tort(TE, G1, G2, N1, N2, alpha, M0, D0):
    return M_ogse_free(TE, G1, N1, TE / N1, M0, alpha * D0) - M_ogse_free(TE, G2, N2, TE / N2, M0, alpha * D0)


def OGSE_contrast_vs_g_mixed(TE, G1, G2, N1, N2, tc, alpha, M0, D0):
    return M_ogse_mixed(TE, G1, N1, TE / N1, tc, alpha, M0, D0) - M_ogse_mixed(TE, G2, N2, TE / N2, tc, alpha, M0, D0)


# ---------------------------------------------------------------------------
# PGSE (mono-exponential) model
# ---------------------------------------------------------------------------

def PGSE_vs_bvalue_exp(bvalue, M0, D0):
    return M0 * np.exp(-bvalue * D0)


