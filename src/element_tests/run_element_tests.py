"""
Validation 1: Single-element NorSand tests.

Runs drained CD triaxial + plane strain for all 8 ψ_0 values,
plus undrained triaxial for ψ_0 = +0.05 and −0.10.

Establishes ε_ref from the triaxial test on ψ_0 = −0.25 (peak q strain).
Extracts φ_mob,TX and φ_mob,PS at ε_ref for every ψ_0.

Outputs:
  outputs/element_tests/phi_mob_table.csv
  outputs/element_tests/element_test_curves.png
  outputs/element_tests/undrained_paths.png
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from src.constitutive.norsand import NorSandParams, NorSandPoint, e_csl

# -----------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------
PARAMS = NorSandParams()
P0 = 32.9          # kPa
MAX_EPS_A = 0.20   # integrate to 20% axial strain
N_STEPS = 500      # strain steps per test

PSI0_VALUES = [-0.25, -0.20, -0.15, -0.10, -0.05, 0.00, +0.05, +0.10]
PSI0_LABELS = ['s0a', 's0b', 's1', 's2', 's3', 's4', 's5', 's6']

# Spec §8 table — void ratios at p'_0 = 32.9 kPa
E0_SPEC = {
    -0.25: 0.566,
    -0.20: 0.616,
    -0.15: 0.666,
    -0.10: 0.716,
    -0.05: 0.766,
     0.00: 0.816,
    +0.05: 0.866,
    +0.10: 0.916,
}


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def run_test(psi0: float, test_type: str) -> NorSandPoint:
    state = NorSandPoint(P0, psi0, PARAMS, test_type=test_type)
    deps = MAX_EPS_A / N_STEPS
    step_fn = {
        'triaxial':            state.step_drained_triaxial,
        'plane_strain':        state.step_drained_plane_strain,
        'undrained_triaxial':  state.step_undrained_triaxial,
    }[test_type]
    for _ in range(N_STEPS):
        step_fn(deps)
    return state


def eps_ref_from_state(state: NorSandPoint):
    """Return (ε_ref, idx) from the triaxial test on the densest sample."""
    q_vals = [h['q'] for h in state.history]
    idx = int(np.argmax(q_vals))
    return state.history[idx]['eps_a'], idx


def interp_phi_at(state: NorSandPoint, eps_target: float) -> float:
    eps_arr = np.array([h['eps_a'] for h in state.history])
    phi_arr = np.array([h['phi_mob'] for h in state.history])
    if eps_target > eps_arr[-1]:
        return float(phi_arr[-1])
    return float(np.interp(eps_target, eps_arr, phi_arr))


def peak_phi(state: NorSandPoint) -> float:
    return max(h['phi_mob'] for h in state.history)


def peak_q(state: NorSandPoint) -> float:
    return max(h['q'] for h in state.history)


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    os.makedirs('outputs/element_tests', exist_ok=True)

    print("Drained triaxial tests...")
    tx_states = []
    for psi0, label in zip(PSI0_VALUES, PSI0_LABELS):
        s = run_test(psi0, 'triaxial')
        tx_states.append(s)
        phi_pk = peak_phi(s)
        q_pk   = peak_q(s)
        print(f"  ψ={psi0:+.2f}  φ_peak={phi_pk:.1f}°  q_peak={q_pk:.1f} kPa"
              f"  e_0={s.history[0]['e']:.4f}")

    print("\nEstablishing ε_ref from ψ_0 = −0.25 triaxial...")
    eps_ref, idx_peak = eps_ref_from_state(tx_states[0])
    print(f"  ε_ref = {eps_ref*100:.3f}%  (step {idx_peak}/{N_STEPS})")

    print("\nDrained plane-strain tests...")
    ps_states = []
    for psi0, label in zip(PSI0_VALUES, PSI0_LABELS):
        s = run_test(psi0, 'plane_strain')
        ps_states.append(s)
        phi_pk = peak_phi(s)
        print(f"  ψ={psi0:+.2f}  φ_peak={phi_pk:.1f}°")

    print("\nUndrained triaxial tests...")
    ud_states = {}
    for psi0 in [+0.05, -0.10]:
        s = run_test(psi0, 'undrained_triaxial')
        ud_states[psi0] = s
        phi_pk = peak_phi(s)
        print(f"  ψ={psi0:+.2f}  φ_peak_eff={phi_pk:.1f}°")

    # -----------------------------------------------------------------------
    # φ_mob table
    # -----------------------------------------------------------------------
    rows = []
    for i, (psi0, label) in enumerate(zip(PSI0_VALUES, PSI0_LABELS)):
        e0 = tx_states[i].history[0]['e']
        phi_tx  = interp_phi_at(tx_states[i], eps_ref)
        phi_ps  = interp_phi_at(ps_states[i], eps_ref)
        phi_pk_tx = peak_phi(tx_states[i])
        phi_pk_ps = peak_phi(ps_states[i])
        rows.append({
            'label':       label,
            'psi_0':       psi0,
            'e_0':         round(e0, 4),
            'eps_ref':     round(eps_ref, 5),
            'phi_mob_TX':  round(phi_tx, 2),
            'phi_mob_PS':  round(phi_ps, 2),
            'phi_peak_TX': round(phi_pk_tx, 2),
            'phi_peak_PS': round(phi_pk_ps, 2),
        })

    df = pd.DataFrame(rows)
    csv_path = 'outputs/element_tests/phi_mob_table.csv'
    df.to_csv(csv_path, index=False)
    print(f"\nφ_mob table → {csv_path}")
    print(df.to_string(index=False))

    # -----------------------------------------------------------------------
    # Spec validation checks
    # -----------------------------------------------------------------------
    print("\n--- Spec §10 checks (triaxial peak) ---")
    expected = {
        -0.25: (46, 48), -0.20: (44, 46), -0.15: (41, 43), -0.10: (35, 37),
        -0.05: (31, 35),  0.00: (29, 33),  0.05: (20, 28),  0.10: (17, 24),
    }
    all_ok = True
    for row in rows:
        psi = row['psi_0']
        phi = row['phi_peak_TX']
        if psi in expected:
            lo, hi = expected[psi]
            ok = lo <= phi <= hi
            status = "OK" if ok else f"WARN (expected {lo}–{hi}°)"
            if not ok:
                all_ok = False
            print(f"  ψ={psi:+.2f}  φ_peak_TX={phi:.1f}°  {status}")
    if all_ok:
        print("  All within spec ranges!")

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    colors = plt.cm.RdYlGn(np.linspace(0.1, 0.9, len(PSI0_VALUES)))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # q vs ε_a (triaxial)
    ax = axes[0, 0]
    for i, (s, lbl, psi0) in enumerate(zip(tx_states, PSI0_LABELS, PSI0_VALUES)):
        eps = [h['eps_a'] * 100 for h in s.history]
        q   = [h['q'] for h in s.history]
        ax.plot(eps, q, color=colors[i], label=f'{lbl} ψ={psi0:+.2f}')
    ax.axvline(eps_ref * 100, color='k', ls='--', lw=1, label=f'ε_ref={eps_ref*100:.2f}%')
    ax.set_xlabel('Axial strain ε_a (%)'); ax.set_ylabel('q (kPa)')
    ax.set_title('Drained Triaxial — q vs ε_a')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # φ_mob vs ε_a (triaxial)
    ax = axes[0, 1]
    for i, (s, lbl, psi0) in enumerate(zip(tx_states, PSI0_LABELS, PSI0_VALUES)):
        eps  = [h['eps_a'] * 100 for h in s.history]
        phi  = [h['phi_mob'] for h in s.history]
        ax.plot(eps, phi, color=colors[i], label=f'{lbl} ψ={psi0:+.2f}')
    ax.axvline(eps_ref * 100, color='k', ls='--', lw=1)
    ax.axhline(PARAMS.phi_c, color='grey', ls=':', lw=1, label=f'φ_c={PARAMS.phi_c}°')
    ax.set_xlabel('Axial strain ε_a (%)'); ax.set_ylabel('φ_mob (°)')
    ax.set_title('Drained Triaxial — φ_mob vs ε_a')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # q vs ε_a (plane strain)
    ax = axes[1, 0]
    for i, (s, lbl, psi0) in enumerate(zip(ps_states, PSI0_LABELS, PSI0_VALUES)):
        eps = [h['eps_a'] * 100 for h in s.history]
        q   = [h['q'] for h in s.history]
        ax.plot(eps, q, color=colors[i], label=f'{lbl} ψ={psi0:+.2f}')
    ax.axvline(eps_ref * 100, color='k', ls='--', lw=1)
    ax.set_xlabel('Axial strain ε_a (%)'); ax.set_ylabel('q (kPa)')
    ax.set_title('Drained Plane Strain — q vs ε_a')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # φ_mob summary scatter
    ax = axes[1, 1]
    phi_tx_list = [interp_phi_at(s, eps_ref) for s in tx_states]
    phi_ps_list = [interp_phi_at(s, eps_ref) for s in ps_states]
    ax.scatter(PSI0_VALUES, phi_tx_list, color='steelblue', marker='o', zorder=3, label='φ_mob,TX')
    ax.scatter(PSI0_VALUES, phi_ps_list, color='tomato',    marker='s', zorder=3, label='φ_mob,PS')
    ax.plot(PSI0_VALUES, phi_tx_list, 'b--', alpha=0.5)
    ax.plot(PSI0_VALUES, phi_ps_list, 'r--', alpha=0.5)
    ax.axhline(PARAMS.phi_c, color='grey', ls=':', label=f'φ_c={PARAMS.phi_c}°')
    ax.set_xlabel('ψ_0'); ax.set_ylabel('φ_mob at ε_ref (°)')
    ax.set_title('φ_mob,TX vs φ_mob,PS at ε_ref')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig('outputs/element_tests/element_test_curves.png', dpi=150)
    plt.close(fig)
    print("\n  → outputs/element_tests/element_test_curves.png")

    # Undrained stress paths
    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 5))
    for psi0, color in zip([+0.05, -0.10], ['tomato', 'steelblue']):
        s = ud_states[psi0]
        p_path  = [h['p'] for h in s.history]
        q_path  = [h['q'] for h in s.history]
        eps_arr = [h['eps_a'] * 100 for h in s.history]
        lbl = f'ψ_0={psi0:+.2f}'
        ax2[0].plot(p_path, q_path, color=color, label=lbl)
        ax2[0].plot(p_path[0], q_path[0], 'o', color=color)
        ax2[1].plot(eps_arr, q_path, color=color, label=lbl)

    p_range = np.linspace(0.1, P0 * 2, 100)
    ax2[0].plot(p_range, PARAMS.M_tc * p_range, 'k--', label=f'M_tc·p')
    ax2[0].set_xlabel("p' (kPa)"); ax2[0].set_ylabel("q (kPa)")
    ax2[0].set_title("Undrained — Effective Stress Path"); ax2[0].legend(); ax2[0].grid(True, alpha=0.3)
    ax2[1].set_xlabel('ε_a (%)'); ax2[1].set_ylabel('q (kPa)')
    ax2[1].set_title("Undrained — q vs ε_a"); ax2[1].legend(); ax2[1].grid(True, alpha=0.3)

    plt.tight_layout()
    fig2.savefig('outputs/element_tests/undrained_paths.png', dpi=150)
    plt.close(fig2)
    print("  → outputs/element_tests/undrained_paths.png")

    print("\nDone.")
    return df, eps_ref


if __name__ == '__main__':
    main()
