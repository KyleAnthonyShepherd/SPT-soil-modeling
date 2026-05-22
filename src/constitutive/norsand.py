"""
NorSand constitutive model (Jefferies 1993; Jefferies & Been 2016).

Sign convention: compression positive for stresses and strains.
Stresses in kPa.

Key equations
-------------
Yield surface : f = η − M_i · [1 − ln(p/p_i)/N] = 0
Image ratio   : M_i = M_tc − χ_tc · ψ_i     (linear state-dilatancy)
State at image: ψ_i = e − e_csl(p_i)
Flow rule     : dε_v^p = D·Λ,  dε_s^p = Λ,  D = M_tc − η  (Rowe)
Hardening     : dp_i/dε_s^p = H · p_i · (M_i − M_tc)
H             : H = H_0 + H_ψ · ψ_i   (decreases toward CSL)

Initialization (isotropic, INSIDE yield surface):
  p_i = p_0    (image cap at current consolidation stress)
  f_0 = −M_i < 0  (elastic interior)
  Soil first yields when stress path hits the yield surface.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

class NorSandParams:
    def __init__(self, phi_c=31.0, Gamma=0.910, lambda_e=0.027,
                 chi_tc=3.5, N=0.30, H_0=100.0, H_psi=200.0,
                 G_ref=50000.0, p_ref=100.0, n_G=0.50, nu=0.20):
        self.phi_c    = phi_c
        self.Gamma    = Gamma
        self.lambda_e = lambda_e
        self.chi_tc   = chi_tc
        self.N        = N
        self.H_0      = H_0
        self.H_psi    = H_psi
        self.G_ref    = G_ref   # kPa
        self.p_ref    = p_ref   # kPa
        self.n_G      = n_G
        self.nu       = nu

        sin_phi       = np.sin(np.radians(phi_c))
        self.M_tc     = 6.0 * sin_phi / (3.0 - sin_phi)
        self.M_te     = 6.0 * sin_phi / (3.0 + sin_phi)


ERKSAK = NorSandParams()


# ---------------------------------------------------------------------------
# Scalar utilities
# ---------------------------------------------------------------------------

def e_csl(p, pr: NorSandParams) -> float:
    return pr.Gamma - pr.lambda_e * np.log(max(p, 1e-9))


def psi_at(e, p, pr: NorSandParams) -> float:
    return e - e_csl(p, pr)


def G_el(p, pr: NorSandParams) -> float:
    return pr.G_ref * (max(p, 1e-9) / pr.p_ref) ** pr.n_G


def K_el(p, pr: NorSandParams) -> float:
    G = G_el(p, pr)
    return 2.0 * G * (1.0 + pr.nu) / (3.0 * (1.0 - 2.0 * pr.nu))


def M_image(e, p_i, pr: NorSandParams) -> float:
    """Image stress ratio M_i = M_tc − χ · ψ_i."""
    psi_i = psi_at(e, max(p_i, 1e-9), pr)
    return pr.M_tc - pr.chi_tc * psi_i


def f_yield(p, q, p_i, M_i, N) -> float:
    """f < 0 → elastic; f = 0 → on surface."""
    if p <= 1e-12 or p_i <= 1e-12:
        return 1e10
    eta = q / p
    return eta - M_i * (1.0 - np.log(p / p_i) / N)


def phi_mob_tx(p, q) -> float:
    """φ_mob (degrees) from triaxial stress invariants."""
    if p < 1e-12:
        return 0.0
    eta = max(q / p, 0.0)
    sin_phi = np.clip(3.0 * eta / (6.0 + eta), -1.0, 1.0)
    return np.degrees(np.arcsin(sin_phi))


# ---------------------------------------------------------------------------
# Triaxial flow-direction helpers
# ---------------------------------------------------------------------------

def flow_m1(D: float) -> float:
    """dε_1^p / dε_s^p  for axisymmetric triaxial."""
    return (3.0 + D) / 3.0


def flow_m3(D: float) -> float:
    """dε_3^p / dε_s^p  for axisymmetric triaxial."""
    return (2.0 * D - 3.0) / 6.0


def drained_tx_modulus(G, K) -> float:
    """E_TX = dσ_1/dε_1 at constant σ_3."""
    lam = K - 2.0 * G / 3.0
    C11 = K + 4.0 * G / 3.0
    C12 = lam                          # = K − 2G/3
    Q   = 2.0 * K + 2.0 * G / 3.0    # C12 + C11
    return C11 - 2.0 * C12 ** 2 / Q


# ---------------------------------------------------------------------------
# Single material-point driver
# ---------------------------------------------------------------------------

class NorSandPoint:
    """
    Single NorSand material point.

    Stress state: (σ_1, σ_3)  compression-positive, σ_2 = σ_3 (triaxial).
    Internal:     e (void ratio), p_i (image mean stress).

    Initialization: isotropic at p0, q0 = 0.
      p_i = p0  →  initial state INSIDE yield surface (elastic interior).
    """

    def __init__(self, p0: float, psi0: float,
                 params: NorSandParams = None,
                 test_type: str = 'triaxial'):
        pr = params or NorSandParams()
        self.pr        = pr
        self.test_type = test_type

        self.sigma1 = float(p0)   # axial effective stress
        self.sigma3 = float(p0)   # radial (= lateral) effective stress
        self.sigma2 = float(p0)   # for plane-strain σ_2 tracking

        self.e   = e_csl(p0, pr) + psi0
        self.p_i = float(p0)       # image cap starts at current p

        # Strain accumulators
        self.eps1  = 0.0
        self.eps3  = 0.0
        self.eps_v = 0.0   # volumetric (compression +)
        self.eps_s = 0.0   # deviatoric shear

        self.history = []
        self._record()

    # ------------------------------------------------------------------
    @property
    def p(self):
        return (self.sigma1 + 2.0 * self.sigma3) / 3.0

    @property
    def q(self):
        return max(self.sigma1 - self.sigma3, 0.0)

    @property
    def eta(self):
        return self.q / self.p if self.p > 1e-12 else 0.0

    # ------------------------------------------------------------------
    def _M_i(self):
        return M_image(self.e, self.p_i, self.pr)

    def _H(self):
        psi_i = psi_at(self.e, self.p_i, self.pr)
        return self.pr.H_0 + self.pr.H_psi * psi_i

    def _f(self):
        return f_yield(self.p, self.q, self.p_i, self._M_i(), self.pr.N)

    # ------------------------------------------------------------------
    def _single_sub_step_tx(self, deps1: float):
        """
        One sub-step of drained triaxial compression (σ_3 = const).

        deps1 : axial strain increment (compression positive).
        """
        pr   = self.pr
        p_c  = self.p     # current mean stress for moduli
        G    = G_el(p_c, pr)
        K    = K_el(p_c, pr)
        lam  = K - 2.0 * G / 3.0
        C11  = K + 4.0 * G / 3.0
        C12  = lam
        Q    = C11 + C12    # = 2K + 2G/3

        # Elastic trial (σ_3 = const → dε_3^e = -C12/Q * deps1)
        deps3_e = -C12 / Q * deps1
        E_TX    = C11 - 2.0 * C12 ** 2 / Q   # drained TX modulus
        dsig1_e = E_TX * deps1
        sig1_tr = self.sigma1 + dsig1_e
        sig3_tr = self.sigma3                  # unchanged

        p_tr  = (sig1_tr + 2.0 * sig3_tr) / 3.0
        q_tr  = max(sig1_tr - sig3_tr, 0.0)
        eta_tr = q_tr / p_tr if p_tr > 1e-12 else 0.0

        # Current M_i (evaluated at pre-step state)
        M_i0   = self._M_i()
        f_tr   = f_yield(p_tr, q_tr, self.p_i, M_i0, pr.N)

        # ----- elastic step -----
        if f_tr <= 1e-10:
            self.sigma1 = sig1_tr
            deps3 = deps3_e
            deps_v = deps1 + 2.0 * deps3
            self.e    -= (1.0 + self.e) * deps_v
            self.eps1  += deps1
            self.eps3  += deps3
            self.eps_v += deps_v
            de_dev = deps1 - deps3
            self.eps_s += abs(de_dev) * 2.0 / 3.0
            return

        # ----- elastoplastic step -----
        # Flow direction in triaxial (Rowe)
        D_tr = pr.M_tc - eta_tr          # dilatancy at trial stress
        m1   = flow_m1(D_tr)             # dε_1^p / dε_s^p
        m3   = flow_m3(D_tr)             # dε_3^p / dε_s^p

        # Additional lateral strain to maintain σ_3=const during plastic flow:
        # dε_3(Λ) = dε_3^e + Λ*(m3 + C12*m1/Q)
        lat_coeff = m3 + C12 * m1 / Q

        # Volumetric strain change per unit Λ (beyond elastic predictor)
        # d(dε_v)/dΛ = 2*lat_coeff + (d/dΛ of axial adjustment)
        # Axial is fixed at deps1, so:
        # dε_v(Λ) = (deps1 + 2*dε_3^e) + 2*lat_coeff*Λ  = dε_v_e + 2*lat_coeff*Λ
        depsv_e = deps1 + 2.0 * deps3_e

        # Hardening constants (use pre-step values, updated after)
        H0     = self._H()
        dpi_dL = H0 * self.p_i * (M_i0 - pr.M_tc)   # dp_i / dΛ

        # Newton solve: f(Λ) = 0
        # σ_1(Λ) = sig1_tr - E_TX*m1*Λ   (constraint-consistent correction)
        # σ_3(Λ) = sig3_tr
        # p_i(Λ) = p_i + dpi_dL*Λ
        # e(Λ)   = e - (1+e)*(depsv_e + 2*lat_coeff*Λ)
        # M_i(Λ) = M_tc - chi*(e(Λ) - e_csl(p_i(Λ)))

        e_now  = self.e
        pi_now = self.p_i

        def state_at(lam):
            s1 = sig1_tr - E_TX * m1 * lam
            s3 = sig3_tr
            if s1 < s3:
                s1 = s3   # clamp (shouldn't happen for small lam)
            pn  = (s1 + 2.0 * s3) / 3.0
            qn  = s1 - s3
            pin = pi_now + dpi_dL * lam
            pin = max(pin, 1e-3)
            depsv = depsv_e + 2.0 * lat_coeff * lam
            en  = e_now - (1.0 + e_now) * depsv
            Mi  = M_image(en, pin, pr)
            return pn, qn, pin, en, Mi

        def residual(lam):
            pn, qn, pin, en, Mi = state_at(lam)
            return f_yield(pn, qn, pin, Mi, pr.N)

        # Newton iterations
        lam = 0.0
        for _it in range(40):
            r  = residual(lam)
            dr = (residual(lam + 1e-9) - r) / 1e-9
            if abs(dr) < 1e-30:
                break
            lam_new = lam - r / dr
            lam = max(0.0, lam_new)
            if abs(r) < 1e-12:
                break

        # Accept corrected state
        pn, qn, pin, en, Mi = state_at(lam)
        self.sigma1 = sig3_tr + qn
        self.sigma3 = sig3_tr
        self.p_i    = pin
        self.e      = en

        deps3   = deps3_e + lat_coeff * lam
        deps_v  = deps1 + 2.0 * deps3
        self.eps1  += deps1
        self.eps3  += deps3
        self.eps_v += deps_v
        de_dev = deps1 - deps3
        self.eps_s += abs(de_dev) * 2.0 / 3.0

    # ------------------------------------------------------------------
    def _single_sub_step_ps(self, deps1: float):
        """
        One sub-step of drained plane-strain compression.
        ε_2 = 0 fixed; σ_3 = 0 (lateral free); σ_2 emerges from constraint.

        Tracks full principal-stress vector [σ_1, σ_2, σ_3].
        """
        pr = self.pr
        G  = G_el(self.p, pr)
        K  = K_el(self.p, pr)
        lam_e = K - 2.0 * G / 3.0
        C = np.array([
            [K + 4*G/3, lam_e,       lam_e      ],
            [lam_e,     K + 4*G/3,   lam_e      ],
            [lam_e,     lam_e,       K + 4*G/3  ]
        ])

        # Constraints: ε_2 = 0 (dε_2 = 0), σ_3 = 0 (dσ_3 = 0)
        # → solve for dε_3 from dσ_3 = 0:
        #   C[2,0]*deps1 + C[2,1]*0 + C[2,2]*deps3 = 0 → deps3 = -C[2,0]/C[2,2]*deps1
        deps3_e = -C[2, 0] / C[2, 2] * deps1
        deps2_e = 0.0

        dsig = C @ np.array([deps1, deps2_e, deps3_e])
        sig_tr = np.array([self.sigma1, self.sigma2, self.sigma3]) + dsig

        p_tr = sig_tr.mean()
        s_tr = sig_tr - p_tr
        q_tr = np.sqrt(1.5 * np.dot(s_tr, s_tr))

        # M_i and yield
        M_i0 = self._M_i()
        f_tr = f_yield(p_tr, q_tr, self.p_i, M_i0, pr.N)

        if f_tr <= 1e-10:
            self.sigma1, self.sigma2, self.sigma3 = sig_tr
            deps_v = deps1 + deps2_e + deps3_e
            self.e    -= (1.0 + self.e) * deps_v
            self.eps1  += deps1
            self.eps3  += deps3_e
            self.eps_v += deps_v
            de_dev = deps1 - deps3_e  # approx
            self.eps_s += abs(de_dev) * 2.0 / 3.0
            return

        # Elastoplastic: radial return in p-q space (plane strain)
        eta_tr = q_tr / p_tr if p_tr > 1e-12 else 0.0
        D_tr   = pr.M_tc - eta_tr
        H0     = self._H()
        dpi_dL = H0 * self.p_i * (M_i0 - pr.M_tc)

        e_now  = self.e
        pi_now = self.p_i

        def residual_ps(lam):
            # Approximate: radial return in p-q space (not constraint-exact)
            p_c = p_tr - K * D_tr * lam
            q_c = q_tr - 3.0 * G * lam
            if q_c < 0:
                q_c = 0.0
            pin = pi_now + dpi_dL * lam
            pin = max(pin, 1e-3)
            depsv = (deps1 + deps3_e) - D_tr * lam   # approx
            en  = e_now - (1.0 + e_now) * (deps1 + deps2_e + deps3_e)
            # (void ratio updated with total elastic vol strain; plastic handled via M_i)
            Mi  = M_image(en, pin, pr)
            return f_yield(p_c, q_c, pin, Mi, pr.N)

        lam = 0.0
        for _it in range(40):
            r  = residual_ps(lam)
            dr = (residual_ps(lam + 1e-9) - r) / 1e-9
            if abs(dr) < 1e-30:
                break
            lam = max(0.0, lam - r / dr)
            if abs(r) < 1e-12:
                break

        p_c = p_tr - K * D_tr * lam
        q_c = max(q_tr - 3.0 * G * lam, 0.0)
        pin = pi_now + dpi_dL * lam
        pin = max(pin, 1e-3)
        deps_v_p = D_tr * lam

        # Reconstruct stress tensor from corrected p, q (preserve direction)
        if q_tr > 1e-12:
            s_c = s_tr * (q_c / q_tr)
        else:
            s_c = np.zeros(3)
        sig_c = s_c + p_c

        self.sigma1, self.sigma2, self.sigma3 = sig_c
        self.p_i = pin
        deps_v = deps1 + deps2_e + deps3_e + deps_v_p
        self.e -= (1.0 + self.e) * (deps1 + deps2_e + deps3_e)  # elastic void ratio
        self.eps1  += deps1
        self.eps3  += deps3_e
        self.eps_v += deps_v
        de_dev = deps1 - deps3_e
        self.eps_s += abs(de_dev) * 2.0 / 3.0

    # ------------------------------------------------------------------
    def _single_sub_step_ud(self, deps1: float):
        """
        One sub-step of undrained triaxial (dε_v = 0).
        σ_3 = const (cell pressure). Excess PWP tracked implicitly via
        effective stress path.
        """
        pr  = self.pr
        G   = G_el(self.p, pr)
        K   = K_el(self.p, pr)

        # Undrained: dε_3 = -deps1/2 (incompressible)
        deps3 = -0.5 * deps1
        deps_s = (2.0 / 3.0) * (deps1 - deps3)   # deviatoric strain = deps1

        # Elastic trial (p unchanged for undrained if purely deviatoric)
        dq_e  = 3.0 * G * deps_s
        sig1_tr = self.sigma1 + dq_e    # dσ_3 effectively 0 in undrained
        sig3_tr = self.sigma3

        p_tr  = (sig1_tr + 2.0 * sig3_tr) / 3.0
        q_tr  = max(sig1_tr - sig3_tr, 0.0)
        eta_tr = q_tr / p_tr if p_tr > 1e-12 else 0.0

        M_i0   = self._M_i()
        f_tr   = f_yield(p_tr, q_tr, self.p_i, M_i0, pr.N)

        if f_tr <= 1e-10:
            self.sigma1 = sig1_tr
            self.eps1  += deps1
            self.eps3  += deps3
            self.eps_v += 0.0   # undrained
            self.eps_s += deps_s
            return

        # Elastoplastic undrained:
        # Undrained constraint: dε_v^total = 0 → dε_v^e = -dε_v^p = -D*Λ
        # → dp^e = K*dε_v^e = -K*D*Λ → p changes despite undrained
        D_tr   = pr.M_tc - eta_tr
        H0     = self._H()
        dpi_dL = H0 * self.p_i * (M_i0 - pr.M_tc)

        e_now  = self.e
        pi_now = self.p_i

        def residual_ud(lam):
            p_c  = p_tr - K * D_tr * lam   # undrained p-shift
            q_c  = max(q_tr - 3.0 * G * lam, 0.0)
            pin  = max(pi_now + dpi_dL * lam, 1e-3)
            Mi   = M_image(e_now, pin, pr)   # e unchanged (undrained)
            return f_yield(p_c, q_c, pin, Mi, pr.N)

        lam = 0.0
        for _it in range(40):
            r  = residual_ud(lam)
            dr = (residual_ud(lam + 1e-9) - r) / 1e-9
            if abs(dr) < 1e-30:
                break
            lam = max(0.0, lam - r / dr)
            if abs(r) < 1e-12:
                break

        p_c  = p_tr - K * D_tr * lam
        q_c  = max(q_tr - 3.0 * G * lam, 0.0)
        pin  = max(pi_now + dpi_dL * lam, 1e-3)

        self.sigma1 = p_c + 2.0 * q_c / 3.0
        self.sigma3 = p_c - q_c / 3.0
        self.p_i    = pin
        # e unchanged (undrained)

        self.eps1  += deps1
        self.eps3  += deps3
        self.eps_v += 0.0
        self.eps_s += deps_s

    # ------------------------------------------------------------------
    def step_drained_triaxial(self, deps1: float, n_sub: int = 10):
        for _ in range(n_sub):
            self._single_sub_step_tx(deps1 / n_sub)
        self._record()

    def step_drained_plane_strain(self, deps1: float, n_sub: int = 10):
        for _ in range(n_sub):
            self._single_sub_step_ps(deps1 / n_sub)
        self._record()

    def step_undrained_triaxial(self, deps1: float, n_sub: int = 10):
        for _ in range(n_sub):
            self._single_sub_step_ud(deps1 / n_sub)
        self._record()

    # ------------------------------------------------------------------
    def phi_mobilised(self) -> float:
        return phi_mob_tx(self.p, self.q)

    def _record(self):
        psi_cur = psi_at(self.e, self.p, self.pr)
        psi_i   = psi_at(self.e, self.p_i, self.pr)
        self.history.append({
            'eps_a':   self.eps1,
            'eps_v':   self.eps_v,
            'eps_s':   self.eps_s,
            'p':       self.p,
            'q':       self.q,
            'e':       self.e,
            'psi':     psi_cur,
            'psi_i':   psi_i,
            'p_i':     self.p_i,
            'M_i':     self._M_i(),
            'phi_mob': self.phi_mobilised(),
            'sigma1':  self.sigma1,
            'sigma3':  self.sigma3,
        })
