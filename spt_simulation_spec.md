# SPT Simulation Specification: Permeability (Fines Proxy) Effect on N

## 1. Research Objective

Quantify how soil permeability — used here as a proxy for non-plastic fines content — affects the SPT blow count N, while holding intrinsic soil strength constant. The expected result is a curve family showing two branches:

- **Loose, contractive soils**: N decreases as k decreases (positive excess PWP buildup reduces effective stress at the tip)
- **Dense, dilative soils**: N increases as k decreases (suction buildup increases effective stress at the tip; this is the Terzaghi-Peck effect)

The deliverable is a defensible drainage-only correction surface N(k, ψ) for use in interpreting field SPT-φ correlations in soils with fines.

**Scope exclusions**: This study deliberately does NOT vary the critical state friction angle φ_c. All cases use the same φ_c. Any apparent peak strength variation must arise from state (ψ) and drainage effects, not from changing the material parameter. Plastic fines (clays) are out of scope.

## 2. Numerical Scheme

**Primary**: Material Point Method (MPM), axisymmetric, coupled u-p (effective stress + pore pressure) formulation. Recommended code: Anura3D v2022 or later.

**Acceptable substitutes** (in order of preference):
- Abaqus/Explicit Coupled Eulerian-Lagrangian (CEL), axisymmetric
- ALE-FEM with adaptive remeshing (Abaqus/Standard, OpenSees, or similar) — only acceptable if remeshing handles 18 inches of penetration cleanly
- DEM with coupled fluid (CFD-DEM) — substantially more expensive, only if continuum methods fail

**Not acceptable**: Standard Lagrangian FEM without remeshing. Pure undrained or pure drained analysis. Mohr-Coulomb constitutive model.

## 3. Geometry

Axisymmetric model (r-z plane, symmetric about z-axis at r=0).

**Soil domain**:
- Radial extent: 500 mm (≈10 outer sampler radii)
- Vertical extent: 1500 mm
- Reference point: sampler tip starts 300 mm below top of soil domain

**Pre-bored hole** (above sampler at start):
- Radial extent: 75 mm
- From z = 0 (top of soil domain) down to z = 300 mm
- Modeled as void (no soil material points)

**Split-spoon sampler** (standard ASTM D1586):
- Outer radius: 25.4 mm (OD = 50.8 mm = 2 inch)
- Inner radius: 17.5 mm (ID = 35.0 mm = 1.375 inch)
- Total length: 600 mm
- Tip geometry: 16° internal bevel on inner edge, simplified to a 1 mm chamfer for meshing
- Material: linear elastic steel, E = 200 GPa, ν = 0.30, ρ = 7850 kg/m³ (rigid is acceptable for Phase 1)
- The sampler is OPEN at the bottom — soil can flow into the inner cavity. Do not enforce a plugged tip.

**Rod above sampler**:
- Length modeled: 1000 mm
- Outer radius: 21.4 mm (AW rod standard)
- Inner radius: 9.5 mm (hollow rod)
- Material: linear elastic steel (same as sampler)
- Connected to sampler at top of spoon with rigid coupling

**Interface**:
- Coulomb friction, δ = (2/3) × φ_c between steel and soil
- Allow separation (no tensile contact)

## 4. Constitutive Model

**Use NorSand** (Jefferies 1993; Jefferies & Been 2016, "Soil Liquefaction: A Critical State Approach", 2nd ed., Chapter 3).

Rationale: φ_c is an explicit input parameter. The state parameter ψ is explicit, which makes the experimental design (sweeping ψ) directly mappable to input. Fewer parameters than SANISAND (8 vs 15) reduce calibration uncertainty. Open implementations exist in MATLAB (Jefferies) and Python (multiple research groups).

**Acceptable substitute**: SANISAND-04 (Dafalias & Manzari 2004) — use if NorSand UMAT is unavailable in target solver. Parameter mapping in §5.

## 5. Material Parameters

### Baseline clean quartz sand (Erksak-330 calibration)

These parameters are held CONSTANT across all simulations.

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Critical state friction angle | φ_c | 31.0 | degrees |
| CSL altitude (e at p'_ref = 1 kPa) | Γ | 0.910 | — |
| CSL slope in e–ln(p') | λ_e | 0.027 | — |
| State-dilatancy parameter | χ_tc | 3.5 | — |
| Volumetric coupling | N | 0.30 | — |
| Plastic hardening modulus multiplier | H_0 | 100 | — |
| Hardening softening parameter | H_ψ | 200 | — |
| Elastic shear modulus reference | G_ref | 50 | MPa |
| Reference pressure for G | p'_ref | 100 | kPa |
| Elastic exponent | n | 0.50 | — |
| Poisson's ratio | ν | 0.20 | — |
| Solid particle density | ρ_s | 2650 | kg/m³ |

Elastic shear modulus: G(p') = G_ref × (p'/p'_ref)^n

### Variables swept across cases

| Variable | Symbol | Range | Units |
|----------|--------|-------|-------|
| Saturated hydraulic conductivity | k | 10⁻² to 10⁻⁸ (7 log-spaced values) | m/s |
| Initial state parameter | ψ_0 | −0.25 to +0.10 (8 values) | — |

Initial void ratio is computed from ψ_0:
```
e_0 = Γ − λ_e × ln(p'_0) + ψ_0
```
where p'_0 is the in-situ mean effective stress at the sampler tip depth (see §6).

### SANISAND-04 equivalent parameters (if substituting)

| SANISAND param | Value | Notes |
|----------------|-------|-------|
| M_c | 1.25 | M_c = 6 sin(31°)/(3−sin(31°)) |
| c (M_e/M_c) | 0.71 | typical for quartz sand |
| λ_c | 0.027 | matches NorSand λ_e |
| e_0 (CSL) | 0.910 | matches Γ |
| ξ | 0.70 | CSL nonlinearity in p^ξ |
| G_0 | 125 | dimensionless |
| ν | 0.20 | |
| h_0 | 7.05 | |
| c_h | 0.968 | |
| n^b | 1.10 | |
| A_0 | 0.704 | |
| n^d | 3.5 | |

## 6. Initial Conditions

**Stress state** (representative of 5 m below water table in a uniform saturated deposit):

- Total unit weight: γ_sat = 20 kN/m³
- Effective unit weight: γ' = γ_sat − γ_w = 10 kN/m³
- Pore pressure at sampler tip depth: u_0 = 50 kPa (assuming WT at original ground surface, 5 m above tip)
- Effective vertical stress at tip: σ'_v0 = 50 kPa
- K_0 = 1 − sin(φ_c) = 0.485
- Effective horizontal stress at tip: σ'_h0 = 24.3 kPa
- Mean effective stress at tip: p'_0 = (σ'_v0 + 2σ'_h0)/3 = 32.9 kPa

**Geostatic gradient** within model domain (depth z measured from top of soil domain):

- σ'_v(z) = (5.0 − 0.3 + z) × γ' for z in meters from soil domain top
- u(z) = (5.0 − 0.3 + z) × γ_w
- σ'_h(z) = K_0 × σ'_v(z)

**Initialization phase**:
1. Apply gravity load with sampler held fixed at z = 300 mm tip depth
2. Equilibrate until kinetic energy < 0.1% of gravitational potential
3. Confirm in-situ stresses match analytical values within 5%

## 7. SPT Loading

**Energy basis**: N_60 (60% energy ratio).

- Theoretical hammer PE: E = m × g × h = 63.5 kg × 9.81 m/s² × 0.762 m = 474.7 J
- Energy delivered (N_60): E_60 = 0.60 × 474.7 = 284.8 J

**Loading approach (simplified, Phase 1)**:

Apply an instantaneous velocity at the top end of the modeled rod, equivalent to a downward-moving rigid mass with KE = 284.8 J.

Effective impact mass = hammer mass + half of modeled rod mass:
- Hammer mass: 63.5 kg
- Rod mass (1 m of AW rod): ρ_steel × π × (r_o² − r_i²) × L = 7850 × π × (0.0214² − 0.0095²) × 1.0 = 9.05 kg
- Effective mass: 63.5 + 9.05/2 ≈ 68.0 kg
- Impact velocity: v_imp = √(2 × 284.8 / 68.0) = **2.89 m/s** downward

Apply v_imp at all material points of the rod top (1 cm slice) for one time step, then release. Let elastic wave propagation handle energy transfer through the rod. Do not constrain rod motion thereafter except by soil-sampler interaction.

**Blow cycling**:
1. Apply blow impulse.
2. Integrate dynamics until kinetic energy in sampler+rod drops below 1% of E_60 (peak settling).
3. Continue simulating consolidation/dynamics for a further 1.5 seconds (representing operator reset time between blows). Pore pressures will partially dissipate at this time scale for high-k cases; for low-k cases they accumulate, which is physical.
4. Record current cumulative tip displacement.
5. If tip displacement < 458 mm AND blow count < 100, apply next blow.

**SPT N convention**:
- N = number of blows required for tip to advance from 152 mm to 458 mm (i.e., the second and third 6-inch increments of an 18-inch drive)
- The first 152 mm (6 inches) is "seating drive" and is NOT counted, per ASTM D1586
- If 100 blows are reached before 458 mm: record as "refusal" with the achieved penetration

## 8. Test Matrix

7 permeability values × 8 state parameters = **56 simulation cases**.

### Permeability sweep

| Case suffix | k (m/s) | Representative soil |
|-------------|---------|----------------------|
| k1 | 1 × 10⁻² | Coarse clean sand / fine gravel |
| k2 | 1 × 10⁻³ | Coarse sand |
| k3 | 1 × 10⁻⁴ | Medium sand |
| k4 | 1 × 10⁻⁵ | Fine sand, trace fines |
| k5 | 1 × 10⁻⁶ | Silty sand, ~15% fines |
| k6 | 1 × 10⁻⁷ | Very silty sand, ~30% fines |
| k7 | 1 × 10⁻⁸ | Sandy silt / silt, 40%+ fines |

### State parameter sweep

| Case suffix | ψ_0 | e_0 (at p'_0 = 32.9 kPa) | Expected behavior |
|-------------|------|---------------------------|---------------------|
| s0a | −0.25 | 0.566 | Extremely dense / OC glacial till regime |
| s0b | −0.20 | 0.616 | Very dense, OC |
| s1 | −0.15 | 0.666 | Dense, strongly dilative |
| s2 | −0.10 | 0.716 | Moderately dense, dilative |
| s3 | −0.05 | 0.766 | Medium dense, modestly dilative |
| s4 | 0.00 | 0.816 | At CSL, neutral |
| s5 | +0.05 | 0.866 | Medium loose, contractive |
| s6 | +0.10 | 0.916 | Loose, strongly contractive |

Note on the densest cases: e_0 = 0.566 (ψ_0 = −0.25) is approaching the minimum void ratio for typical quartz sand (e_min ≈ 0.55), representing overconsolidated glacial till or heavily compacted fill. Confirm e_0 > e_min for the specific sand parameter set being used; if not, reduce |ψ_0| for the densest case.

### Run priority

If computational budget is limited, run in this order:

**Phase 1A (anchor — 1 run)**: s4-k3 (CSL, medium sand) — verifies the implementation against expected N ≈ 10–15 for clean medium sand at moderate density.

**Phase 1B (branch identification — 14 runs)**:
- All 7 k-values at ψ = −0.10 (s2) — dense branch
- All 7 k-values at ψ = +0.05 (s5) — loose branch
- Plot N vs log(k) for both. The dense series should rise with decreasing k; the loose series should fall.

**Phase 1C (very-dense regime — 14 runs)**:
- All 7 k-values at ψ = −0.20 (s0b) — very dense, expect very high N
- All 7 k-values at ψ = −0.25 (s0a) — glacial till regime, refusal possible for low k
- These are the cases most likely to hit refusal (100 blows); document refusal patterns carefully.

**Phase 1D (fill — 27 runs)**: Remaining cases to populate the full N(k, ψ) surface.

## 9. Outputs

### Per-blow records (time series CSV per case)

| Column | Description | Units |
|--------|-------------|-------|
| blow_n | Blow number | — |
| t_blow | Simulation time at blow | s |
| disp_tip | Cumulative tip displacement | mm |
| dpen | Penetration this blow | mm |
| F_tip_peak | Peak axial force on sampler tip during blow | kN |
| F_shaft_peak | Peak shaft friction force during blow | kN |
| u_excess_tip | Peak excess PWP at r = 30 mm from axis, z = mid-sampler | kPa |
| u_excess_postblow | Excess PWP at same point after 1.5 s settle | kPa |
| e_avg_tip | Volume-averaged void ratio in 50 mm radius sphere centered on tip | — |
| sigma_v_eff_tip | Vertical effective stress at tip after blow | kPa |

### Case summary

For each of the 56 cases:
- N (blows from 152 to 458 mm), or "refusal at X blows / Y mm"
- Time-series CSV (above)
- Final void ratio contour plot
- Tip force vs cumulative displacement plot
- Excess PWP vs blow number plot

### Per-ψ_0 element-test extraction (8 cases, run once per ψ_0)

Two drained element-test simulations per ψ_0 at p'_0 = 32.9 kPa with the in-situ e_0:
- **Triaxial compression (CD)** — axisymmetric, σ'_2 = σ'_3
- **Plane strain compression (CD)** — ε_2 = 0, σ'_2 free

Record the full q vs ε_1 curve and σ'_1, σ'_3 history.

Reference strain ε_ref is defined from the **triaxial** test on the densest case (ψ_0 = −0.25): the axial strain at which q reaches its maximum. This is a single number — pin it down once, then use it for all subsequent extractions across both test types and all ψ_0.

For each ψ_0, evaluate at ε_1 = ε_ref:
- φ_mob,TX = arcsin[(σ'_1 − σ'_3) / (σ'_1 + σ'_3)] from the triaxial test
- φ_mob,PS = arcsin[(σ'_1 − σ'_3) / (σ'_1 + σ'_3)] from the plane strain test

These give two columns of 8 values each (one per ψ_0), constant across k.

### Study deliverables

- Master CSV: rows = 56 cases, columns = (k, ψ_0, e_0, N, refusal_flag, mean_F_tip_per_blow, ...)
- φ_mob CSV: rows = 8 ψ_0 values, columns = (ψ_0, e_0, ε_ref, φ_mob_TX, φ_mob_PS, φ_peak_TX, φ_peak_PS)
- Plot 1: N vs log(k) — one line per ψ_0
- Plot 2a: φ_mob,TX vs N, with Wolff 1989 overlaid (primary comparison — Wolff is triaxial-calibrated)
- Plot 2b: φ_mob,PS vs N, same overlay (for strip-footing applications; expected to lie above Plot 2a)
- Plot 3: 2D contour of N over (log k, ψ_0) space
- Plot 4: Excess PWP at tip at blow 1 vs log(k) — confirms drainage regime

## 10. Validation

Complete BEFORE running the full sweep.

### Validation 1: Single-element response and ε_ref establishment

This is also the source of the φ_mob,TX and φ_mob,PS columns in the deliverables.

Run drained CD compression on a single material point at p'_0 = 32.9 kPa for **all 8 ψ_0 values**, in both triaxial and plane strain geometry. Use the e_0 corresponding to each ψ_0 from the §8 table.

From the **triaxial test on ψ_0 = −0.25**, identify the axial strain at peak q. This is ε_ref — the single reference strain used for all φ_mob extractions across the study.

Expected results (triaxial):
- ψ_0 = −0.25: clear peak, φ_peak ≈ 46–48°, ε_ref likely 2–4%
- ψ_0 = −0.20: clear peak, φ_peak ≈ 44–46°
- ψ_0 = −0.15: peak, φ_peak ≈ 41–43°
- ψ_0 = −0.10: modest peak, φ_peak ≈ 35–37°, softens toward φ_c = 31°
- ψ_0 = −0.05: weak peak ≈ 33°
- ψ_0 = 0.00: monotonic toward φ_c with no clear peak
- ψ_0 = +0.05: strain-hardens, φ_mob at ε_ref well below φ_c (~22–26°)
- ψ_0 = +0.10: strongly contractive, φ_mob at ε_ref ~ 19–22°

Expected results (plane strain): peak values approximately 5–8° higher than triaxial across the range; φ_peak,PS for ψ_0 = −0.25 likely 52–55°.

Also run undrained triaxial at ψ_0 = +0.05 and ψ_0 = −0.10:
- Undrained, ψ_0 = +0.05: clear effective stress path with PWP buildup, may reach quasi-steady state (flow liquefaction)
- Undrained, ψ_0 = −0.10: dilative path, suction develops, deviator stress rises monotonically

If any of these deviate significantly, parameters or implementation are wrong.

### Validation 2: Quasi-static drained CPT analog

Push the sampler at constant rate 1 mm/s through 200 mm of penetration. Force fully drained by setting k = 10⁻¹ m/s.

Expected for ψ_0 = −0.05 at the assumed stress level:
- Steady-state tip resistance q_t = 6–10 MPa
- Tip resistance correlation Q_t = (q_t − σ_v) / σ'_v ≈ 120–200
- This is in the range expected for medium-dense clean sand at σ'_v = 50 kPa

### Validation 3: Cavity expansion check

Compute the spherical cavity limit pressure using Carter, Booker & Yeung (1986) for critical state soil with the same parameters. Use the relation q_t ≈ (1 + I_rr/3) × p_limit (Vesic 1972 shape factor, with rigidity index I_rr from G and undrained-equivalent strength).

Expected agreement with Validation 2 result: within ±30%. Cavity expansion is a known approximation; this is a sanity check, not a tight bound.

### Validation 4: Energy audit

For one blow at any case:
- Sum kinetic energy delivered to sampler+rod system at impact
- Compare to E_60 = 284.8 J
- Should match within 5% (residual is wave dispersion before peak)

## 11. Pitfalls and Implementation Notes

1. **Time step**: Coupled u-p schemes are stable only if Δt ≤ min(h_min/c_p, h_min²/(2c_v)), where c_p is P-wave speed and c_v is consolidation coefficient. For low k, the diffusion limit may force Δt < 10⁻⁶ s. Consider an implicit scheme for the pore pressure DOF if explicit time-stepping becomes impractical.

2. **Numerical damping**: Some artificial damping is needed in MPM to suppress spurious oscillations. Use Rayleigh damping with α_M = 0.05/T_blow, where T_blow ≈ 5 ms is the blow duration. Do not exceed this — excessive damping will artificially reduce N.

3. **Mesh/particle resolution**: At least 4 material points per sampler radius near the tip (i.e., characteristic spacing ≤ 6 mm in the tip region). Coarser away from the tip is fine.

4. **Boundary reflections**: P-wave speed in saturated sand ≈ 1500 m/s (water-saturated). A reflection from the 500 mm radial boundary returns in ~0.67 ms. This is comparable to blow duration. Use absorbing boundaries (viscous dashpots) on the radial and bottom faces. Required, not optional.

5. **Sampler plugging**: Track soil that enters the inner cavity. If it reaches the top of the inner spoon before 458 mm of penetration, the soil column inside the spoon will start to lock up and increase resistance. This is physical and should be allowed. Do not artificially eject soil from the cavity.

6. **Multiple-blow soil evolution**: By blow 20, soil ahead of the sampler will be densified relative to initial conditions. This is part of why N increases with depth and is part of the physical response. Do not reset state between blows.

7. **Sign convention for ψ_0**: Negative ψ = denser than critical. Positive ψ = looser than critical. Confirm this matches the constitutive model implementation; some codes invert the sign.

8. **Refusal handling**: If a case reaches 100 blows without 458 mm penetration, record the case as "refusal" but also report the partial N. Do not extrapolate. For very dense / very low-k cases this may happen and is meaningful.

9. **No fitting to field data**: Do not tune any parameter to match field N-φ correlations. The point of this study is to predict the deviation independently and compare. Tuning would defeat the purpose.

10. **Friction angle extraction (replacing earlier "peak φ" procedure)**: For each ψ_0, run drained CD element tests in triaxial and plane strain geometry at p'_0 = 32.9 kPa with the in-situ e_0. Establish ε_ref from the triaxial test on the densest case (ψ_0 = −0.25) — the axial strain at peak q. Then for every ψ_0, evaluate the mobilized friction angle at ε_1 = ε_ref under both test geometries. This gives φ_mob,TX and φ_mob,PS as functions of ψ_0 only — not of k. The k-driven variation in N at constant φ_mob is the central result of the study.

    Why not back-solve from the SPT itself: bearing capacity and cavity expansion equations assume a Mohr-Coulomb failure surface controlled by the same friction angle they're solving for. Back-solving simply reproduces whichever φ would have generated the observed q_t/σ'_v — it doesn't independently estimate φ. The displacement-based element-test definition is operationally meaningful (it ties φ_mob to a tolerable strain that matches what bearing capacity equations were calibrated against) without that circularity.

    Why ε_ref is set by the densest case: at the displacement that pushes the densest sample to its peak, looser samples are still strain-hardening below their critical state strength. This naturally produces lower φ_mob for loose soils — the regime where bearing capacity practice uses 20–25°. It's the same logic used implicitly in DST interpretation.

    Triaxial is primary (Wolff is triaxial-calibrated). Plane strain is secondary, reported for strip-footing applications where the higher mobilized angle is defensible.

## 12. Reporting

Deliver:

1. The 56-case SPT results table (CSV)
2. The 8-row element-test extraction table (φ_mob,TX and φ_mob,PS per ψ_0)
3. The five study plots (Plot 1, 2a, 2b, 3, 4) as PDF or PNG
4. A short narrative (2 pages max) covering:
   - Did the dense branch show N rising with decreasing k? By how much?
   - Did the loose branch show N falling with decreasing k? By how much?
   - At what ψ does the response flip (transition between branches)?
   - What is the magnitude of departure from clean-sand correlations attributable purely to drainage?
   - How much higher is φ_mob,PS than φ_mob,TX across the ψ range? Is there a stable offset, or is it state-dependent?
   - Did any glacial-till regime cases (ψ_0 = −0.20, −0.25) hit refusal? At what k?
   - Where validation cases failed, what was the discrepancy and why?
5. Validation case outputs (8 × 2 element-test curves with ε_ref marked, quasi-static CPT response, cavity expansion comparison, energy audit)
6. Input files for all 56 SPT cases and 16 element-test cases

## References

- Carter, J.P., Booker, J.R., Yeung, S.K. (1986). "Cavity expansion in cohesive frictional soils." *Géotechnique* 36(3): 349–358.
- Cubrinovski, M., Ishihara, K. (2002). "Maximum and minimum void ratio characteristics of sands." *Soils and Foundations* 42(6): 65–78.
- Dafalias, Y.F., Manzari, M.T. (2004). "Simple plasticity sand model accounting for fabric change effects." *Journal of Engineering Mechanics* 130(6): 622–634.
- Jefferies, M.G. (1993). "NorSand: a simple critical state model for sand." *Géotechnique* 43(1): 91–103.
- Jefferies, M., Been, K. (2016). *Soil Liquefaction: A Critical State Approach*, 2nd ed., CRC Press.
- Vesic, A.S. (1972). "Expansion of cavities in infinite soil mass." *Journal of the Soil Mechanics and Foundations Division* 98(SM3): 265–290.
- Wolff, T.F. (1989). "Pile capacity prediction using parameter functions." *ASCE Geotechnical Special Publication* No. 23: 96–106.
- Yu, H.S., Houlsby, G.T. (1991). "Finite cavity expansion in dilatant soils: loading analysis." *Géotechnique* 41(2): 173–183.
