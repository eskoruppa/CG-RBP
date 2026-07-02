# Simulation Methodology for Torque-Induced DNA Melting as Barrier Crossing in the CG-RBP LAMMPS Model

**Scope.** This document specifies the *simulation methodology* — ensembles, collective variables (CVs), free-energy / rate methods, and a concrete prioritized protocol — to study torque-induced duplex melting as a Kramers barrier-crossing problem in the existing SE(3) coarse-grained rigid-base-pair (CG-RBP) LAMMPS force field at ~10 bp / bead resolution. It assumes the two companion design documents (Challenge 1: unwrapped twist; Challenge 2: bimodal double-well energy) and grounds every recommendation in the master literature review (sections cited as **[LR §n]**).

It is written to be implementation-ready against the actual code: `fix_rbp_lrf.cpp/.h`, `so3.h`, `bond_rbp*.cpp`, `dihedral_rbp.cpp`, `angle_rbp.cpp`, `parse_rbp.h`.

---

## 0. Executive summary / prioritized plan

The single most on-target template for this entire program is **Sicard, Destainville, Manghi 2015** (J. Chem. Phys. 142:034903) [LR §7]: a mesoscopic rotating-strand model with a distance-dependent torsional modulus, well-tempered metadynamics with **twist as the rate-limiting CV**, yielding explicit ~22 kBT opening / ~13 kBT closure barriers for a 10-bp bubble. The methodology below adapts that recipe to the SE(3) CG-RBP bead model.

**Prioritized phases (each gated on the prior validating):**

| Phase | Goal | Method | Validation target |
|---|---|---|---|
| **P0** | Equilibrium duplex elasticity reproduced | unbiased NVE/Langevin, free ends | A≈45–50 nm, C≈95–110 nm [LR §1,§C3] |
| **P1** | Linear torque response (single well) | constant applied torque on terminal triad; measure ⟨Tw⟩ | τ = C_eff·2πσ/L slope [LR §2; Marko 2007] |
| **P2** | 1-D PMF along unwrapped twist (the headline result) | umbrella sampling along Φ (winding-aware CV) + MBAR | three-regime G(σ): elastic / transition / flat [LR §7; Liebl-Zacharias] |
| **P3** | Melting torque & coexistence | Maxwell/lever-rule on G(σ); or constant-torque scan | τ_melt ≈ −10 to −11 pN·nm [LR §1,§C5] |
| **P4** | Nucleation barrier & rate | 2-D metadynamics (twist, melted-step count) + Kramers/MFPT; FFS cross-check | ΔF_op≈22, ΔF_cl≈13 kBT [LR §7] |
| **P5** | Sequence dependence & hysteresis | AT-rich vs GC-rich; torque-ramp at varying rates | AT-rich nucleates first [LR §3,§4]; hysteresis loop area → rate |

**Two methodology decisions that pervade everything:**
1. **Never bias a raw rotation-matrix twist.** PLUMED's `TORSION` and any `acos`-based angle are wrapped to (−π,π] [LR §A, §7 caveat]. All biasing/CV machinery consumes the **unwrapped winding accumulator** from Challenge 1 (exposed as a per-atom property and a global `compute`), or the **periodicity-free discrete melted-step count**.
2. **Ensemble duality.** Run *both* a **constant-applied-torque** ensemble (Gibbs/Legendre, tilt `−τ·Tw`) and a **fixed-linking-number** ensemble (clamped terminal triads, measure conjugate torque). Coexistence (degenerate basins) defines τ_c; the two ensembles must agree on τ_c via the Legendre transform [LR §2; Moroz-Nelson; cgNA+min].

---

## 1. Ensembles

### 1.1 Integrator layer

Beads are `atom_style ellipsoid` rigid bodies (the LRF fix already casts to `AtomVecEllipsoid` and reads `bonus[n].quat`). Two viable integrators:

- **`fix nve/asphere` + `fix langevin` (+ `fix langevin` angular)** — standard LAMMPS ellipsoid integration. Simplest; thermostats translational and rotational DOF. Adequate for free-energy sampling where only equilibrium distributions matter.
- **Quaternion-Langevin geometric integrator** (`fix nve/dotc/langevin` from the CG-DNA/oxDNA package; Davidchack-Ouldridge-Tretyakov 2017; Henrich 2018 [LR §4,§7]) — exactly conserves |q|, OU angular-momentum noise, anisotropic friction, correct canonical SE(3) sampling. **Preferred** because (i) it is the validated engine for rigid-body DNA in LAMMPS and (ii) its torque convention is documented and matches the body-frame torques the RBP styles already write to `atom->torque`.

> **Critical note (carries from Challenge 1):** the integrator advances orientation in SO(3) and **does not store winding**. The unwrapped twist must be reconstructed externally each step by the LRF fix (see §2.1). The `fix rbp/lrf` already runs `PRE_FORCE` every step, computes per-bond `fwd_euler` (= Ω or Φ_d) keyed to the left atom, and forward-comms it to ghosts — it is the natural and *only* correct home for the winding counter.

### 1.2 Geometry, anchoring, and the helical axis

- **Topology:** a linear backbone of N beads (one RBP bond per consecutive pair, left→right 5'→3' orientation — the LRF fix's orientation requirement). Optionally closed (minicircle) using the `closed` metadata flag [parse_rbp.h] for a strictly-defined Lk.
- **Helical axis ẑ:** define by the two clamped ends. Apply tension F along ẑ to suppress writhe so that imposed turns go into twist (the experimental regime above F_c≈0.5 pN where melting, not plectoneme, absorbs undertwist [LR §1; Strick 1998; Matek 2015]).
- **Anchoring options** (terminal triads "clamped"):
  - **Bead 0 fully fixed:** `group anchor id 1`; `fix setforce anchor 0 0 0` + zero its torque (a `fix` that zeroes `torque[i]`), and pin its quaternion (a small custom `fix freeze/quat`, or `fix rigid` on a 1-bead group). This fixes the lab-frame reference triad T_0.
  - **Bead N−1 = the "rotor"/handle:** force along ẑ via `fix addforce`; torque/turns about ẑ via the mechanisms in §1.3–1.4.

### 1.3 Constant-applied-torque ensemble (Gibbs / Legendre)

Add a constant external torque τ·ẑ to the terminal handle bead about the helical axis. This is the Bell/Legendre tilt of the landscape by `−τ·Tw` [LR §2,§6; Bell 1978; Marko 2007].

Two equivalent implementations:

1. **Native LAMMPS:** `fix addtorque` (or `fix addtorque/atom`) on the handle group, vector `0 0 τ`. `fix addtorque` applies a torque to a group's COM; for a single terminal bead use `fix addtorque/atom`. Time-dependence (torque ramps for hysteresis, P5) via an `equal`-style variable.
   ```
   variable tau    equal  v_tau0 + v_taurate*time
   fix      hold   handle addtorque/atom 0.0 0.0 v_tau
   ```
2. **Energy-term implementation (recommended for free-energy consistency):** add a potential `E_τ = −τ·Φ` where Φ is the **unwrapped global twist** CV (§2.1). The conjugate force is a constant torque about ẑ distributed to the chain via the same `Jinvtp`→body-torque map the bond style already uses. This guarantees the bias enters the PMF reweighting cleanly (the applied-torque term is exactly the linear tilt in the energy forms B1/B2 of Challenge 2) and avoids COM/handle subtleties. Implement as a small `fix rbp/torque` that, given the per-atom twist gradient `∂Φ/∂q_i` already available from the LRF Jinvtp, adds `−τ ∂Φ/∂(orientation)` to `atom->torque`.

**Measuring the conjugate variable (twist/turns):** record the unwrapped `Φ = Σ_i twist_i` (§2.1) each frame; ⟨Φ⟩(τ) is the torque–turns response curve (validation against the "hat curve" / linear branch [LR §1,§5]).

### 1.4 Fixed-linking-number ensemble (fixed turns)

Impose ΔLk = n exactly (the canonical magnetic-tweezers ensemble; Strick 1996 [LR §1]) and **measure** the conjugate torque.

- **Clamp the terminal triad to a prescribed rotation angle.** Rotate the handle bead's body frame by 2πn about ẑ relative to bead 0 and hold it: a stiff harmonic angular restraint to a target quaternion, or `fix rigid` on the 1-bead handle group with a prescribed orientation. For a ramp, advance the target angle slowly (quasi-static; Deufel 2007 fixed-Lk by polarization rotation [LR §1]).
- **Constrain the CV instead of the handle (preferred, ensemble-clean):** a stiff harmonic restraint on the **global twist CV** Φ (or on Lk = Tw + Wr, §2.3) to target Φ* = 2π·(turns):
  `E_Lk = ½ k_Lk (Φ − Φ*)²`, k_Lk large.
  This is the LAMMPS realization of the cgNA+min "torque as a Lagrange multiplier conjugate to integrated twist" [LR §7].
- **Conjugate torque readout:** τ = −⟨∂F/∂Φ⟩. In the stiff-restraint method, τ = ⟨k_Lk (Φ* − Φ)⟩ (the mean restraint force on the CV) — exactly the rotor-bead `τ = k_ang·⟨Δφ⟩` readout [LR §1; Bryant 2003]. Sweep Φ* to map τ(Φ*): the plateau is τ_melt.

**Why both ensembles:** in the two-state region the fixed-Lk ensemble shows **backbending / negative torsional compliance** (a coexistence signature [LR §1; Oberstrass 2013]) whereas the fixed-torque ensemble hops between basins. The fixed-Lk τ(Lk) plateau and the fixed-torque P(open)(τ) step must reconcile via the Maxwell construction (§3.3).

### 1.5 Writhe control

To keep the experiment in the **twist-dominated** regime (so undertwisting drives melting rather than plectoneme formation), apply sufficient tension (F ≳ 0.5–2.5 pN [LR §1,§4; Matek 2015]) and monitor writhe Wr each frame (§2.3). If Wr grows (buckling), either raise F or accept plectoneme/bubble competition as a separate study. Marko 2007's coexistence picture assumes the stretched branch [LR §2].

---

## 2. Collective variables

All CVs are computed **on the fly** by extending `fix rbp/lrf` (which already has the per-junction `fwd_euler`, ghost comm, and grow/copy/exchange machinery) and exposed to LAMMPS as `compute`s / per-atom properties for PLUMED or native biasing.

### 2.1 Unwrapped (winding-aware) twist — primary CV

**Definition.** Per step i, the twist is the third component of the Y-convention deformation already computed: `t_i^bare = fwd_euler[id1][2]` (= Φ_d,3, the deformation twist that has the intrinsic ~180° branch removed). The **unwrapped per-step twist** carries an integer winding counter w_i:

```
t_i^unwrapped(frame k) = t_i^bare(k) + 2π · w_i
```

**Unwrapping rule (per step, per frame):**
```
Δ = t_i^bare(k) − t_i^bare(k−1)
if Δ >  π:  w_i -= 1
if Δ < −π:  w_i += 1
t_i^unwrapped = t_i^bare(k) + 2π·w_i
```

**Global twist CV (the reaction coordinate):**
```
Φ = Σ_{i=0}^{N-2} t_i^unwrapped
```

**Implementation in `fix rbp/lrf` (concrete):**
- Add a persistent per-atom integer array `wind` (keyed to the left atom, exactly like `fwd_euler`) and a per-atom previous-twist scalar `tw_prev`.
- Allocate/grow them in `grow_arrays`, copy in `copy_arrays`, zero in `set_arrays`, and **add to `pack_exchange`/`unpack_exchange` and `pack_forward_comm`** so winding survives domain migration and reaches ghosts (the same pattern already used for `fwd_euler`). This is the load-bearing extension — the fix currently has NO `write_restart`/`read_restart`, so **add them** to persist `wind`/`tw_prev` across restarts (otherwise the absolute winding number is lost on restart; see §6 pitfalls).
- In `compute_lrf` Phase 2, immediately after `fwd_euler[id1][2]` is set, run the unwrapping rule against `tw_prev[id1]`, update `wind[id1]`, store the unwrapped value, and set `tw_prev[id1]`.
- Expose Φ via a small `compute rbp/twist` (a `Compute` with `scalar_flag=1`) that sums `t_i^unwrapped` over owned left atoms and `MPI_Allreduce`s. Per-step twist exposed as a per-atom `compute` for the melted-step CV.

```cpp
// added to FixRBPLRF (sketch)
double *tw_prev;   // bare twist at previous frame, per left atom
int    *wind;      // winding integer, per left atom
// in compute_lrf Phase 2, after euler[2] computed for Y convention:
double tb = euler[2];
double d  = tb - tw_prev[id1];
if      (d >  MY_PI) wind[id1] -= 1;
else if (d < -MY_PI) wind[id1] += 1;
tw_unwrapped[id1] = tb + 2.0*MY_PI*wind[id1];
tw_prev[id1] = tb;
```

**Dump-interval caveat.** The unwrapping is unambiguous only if per-step bare twist changes < π between *consecutive PRE_FORCE calls* (i.e. every MD step). Because the fix runs every step (not every dump), this is essentially always satisfied at 10-bp resolution and a sane timestep. Choose Δt so a single step never untwists a junction by > π (Challenge-1 caveat). [LR §A; Brackley 2014; Bergou 2008.]

**Fallback (robustness):** parallel-transport (Bishop-frame) twist (Discrete Elastic Rods; Bergou 2008 [LR §A2]) if per-step unwrapping proves fragile in tightly writhed configs. Not needed at the recommended tensions.

### 2.2 Number of melted steps / bubble size — periodicity-free CV

This CV sidesteps angle periodicity entirely and is the natural FFS/metadynamics order parameter [LR §4,§7; oxDNA-FFS; Matek 2012].

**Definition.** A step i is "melted" if its deformation sits in the soft/melted basin of the Challenge-2 double-well. Two equivalent indicators:
- **Energetic:** the softmax weight of the melted basin from the log-sum-exp energy (B1):
  `m_i = exp(−ρ(E_melt,i+ΔG_i)) / [exp(−ρ E_dup,i) + exp(−ρ(E_melt,i+ΔG_i))] ∈ [0,1]`.
  This is *already computed* inside the double-well bond style — expose it as a per-bond/per-atom scalar.
- **Geometric (cheaper, smooth):** a sigmoid on the unwrapped step twist relative to the crossover twist t*:
  `m_i = ½[1 − tanh((t_i^unwrapped − t*)/δ)]`, with t* halfway between duplex (~+0.6 rad/bp·×bp) and melted (~−2π/10 bp) twist.

**Bubble CVs:**
- `M = Σ_i m_i` — number (or fraction) of melted steps (the global denaturation order parameter; Michieletto 2017 φ-field [LR §4]).
- Largest contiguous run of m_i>½ — bubble size (for nucleation-vs-growth analysis).

Expose `M` as `compute rbp/melted` (sum + Allreduce). Use `M` as the **second metadynamics CV** and as the **FFS order parameter** (integer-valued count is ideal for FFS interfaces λ_i, no periodicity issue).

### 2.3 Lk = Tw + Wr topological check (independent cross-validation)

Compute writhe Wr per frame from bead centerlines and check Lk = Tw + Wr against the imposed Lk each frame [LR §5; White; Fuller; Klenin-Langowski].
- **Tw** = Φ/2π (from §2.1).
- **Wr**: Klenin-Langowski exact polygonal solid-angle double sum (O(N²), robust everywhere) for validation runs; Fuller single-integral relative to ẑ for speed (`dWr = (1/2π)(1−cosθ)dφ`), **switching to the full sum near antipodal points** (t·ẑ→−1, plectonemes) where Fuller fails [LR §5; Neukirch-Starostin 2008].
- This is a *diagnostic / optional CV*, not the primary RC. For closed minicircles Lk is a strict integer invariant and gives a hard correctness check; for open clamped chains use Fuller relative writhe + explicit end framing.

Implement as a standalone `compute rbp/writhe` reading `atom->x` of the backbone (gather to root for the O(N²) sum, or a `fix` that maintains it). At 10-bp/bead with N~few hundred beads the O(N²) cost is negligible relative to force evaluation.

### 2.4 Native compute vs PLUMED

- **Native computes (recommended primary):** `compute rbp/twist` and `compute rbp/melted` give exact, winding-aware values. Bias them with LAMMPS `fix restraint`/`fix spring` (umbrella) or `fix smd`. No periodicity hazard because the CV is already unwrapped.
- **PLUMED (`fix plumed`):** preferred for metadynamics/OPES, MBAR-style reweighting, replica/bias-exchange, and PATHMSD [LR §7; Bonomi 2019]. **Do not use PLUMED's `TORSION`** (wrapped to (−π,π]). Instead:
  - feed PLUMED the native unwrapped CV via a LAMMPS per-atom property and PLUMED `CUSTOM`/`MATHEVAL`/`COMBINE`, **or**
  - bias the **melted-step count M** (periodicity-free) directly with PLUMED collective variables built from per-bead positions/distances. The latter is the clean route and is exactly how oxDNA-FFS biases bp count [LR §7].

---

## 3. Free-energy and rate methods

### 3.1 Method selection (matched to the physics)

| Question | Method | CV(s) | Source |
|---|---|---|---|
| Shape of the twist PMF (the double well) | **Umbrella sampling + MBAR** | Φ (unwrapped) | [LR §7; Torrie-Valleau; Shirts-Chodera] |
| 2-D landscape, nucleation barrier | **Well-tempered metadynamics / OPES** | (Φ, M) | [LR §7; Laio-Parrinello; Sicard 2015] |
| Melting rate (no good RC assumed) | **Forward Flux Sampling** | M (interfaces λ_i = melted-step count) | [LR §7; Allen 2005; oxDNA-FFS] |
| Is Φ a good RC? | **Committor / TPS spot-check** | — | [LR §7; Dellago 1998] |
| Optimal path + D(α) for 1-D rate | **String method in CVs** | (Φ, M) plane | [LR §7; Maragliano 2006] |
| Assemble rate from short runs | **MSM / PCCA+** | featurized (Φ, M) | [LR §7; Chodera-Noé] |

### 3.2 Umbrella sampling along twist (Phase P2 — the headline calculation)

This directly mirrors the all-atom blueprint of Liebl-Zacharias (RC = unwinding angle/σ, umbrella + WHAM) [LR §7] and Wereszczynski-Andricioaei (umbrella under applied torque) [LR §2], at 10-bp/bead.

- **Windows:** stratify Φ from duplex (Φ≈0) past one full untwist turn (Φ ≈ −2π·N_bubble/N) into the melted basin. Spacing ≈ 0.2–0.5 turn; **denser near the barrier** (place windows by an initial coarse scan / WHAM error feedback [LR §7; Kumar 1992]).
- **Bias:** harmonic restraint `½k_u(Φ−Φ_0)²` on the native unwrapped CV, k_u chosen so adjacent windows overlap (target ~20–30% histogram overlap).
- **Reweighting:** **MBAR** (binless, multidimensional-ready, rigorous error bars on barrier height [LR §7; Shirts-Chodera 2008]) — preferred over WHAM because P2→P4 extends to 2-D (Φ,M).
- **Output:** G(Φ) (equivalently G(σ), σ = ΔLk/Lk0). Expect the **three-regime shape** [LR §7; Liebl-Zacharias]: (I) harmonic duplex elasticity → (II) abrupt transition/barrier → (III) flat, low-slope melted branch. The flat region's slope is the melting torque (§3.3).

### 3.3 Extracting the melting torque and coexistence

The applied/fixed-torque duality and Maxwell construction [LR §2,§6; Marko 2007; Marko-Neukirch 2013]:

- **From the PMF:** τ(Φ) = dG/dΦ (in units where Φ is the integrated twist; convert per-turn). The **coexistence torque τ_c** is the constant slope of the common tangent (Maxwell construction) between the duplex and melted basins of G(Φ). Equivalently, τ_c is the τ at which the tilted PMF `G(Φ) − τ_c Φ` has two degenerate minima.
- **Calibration constraint (must hold):**
  `τ_c = (ΔG_melt − ΔG_dup)/(Tw0_dup − Tw0_melt)`  [LR §6 B6; Marko-Neukirch]
  must equal **−10 to −11 pN·nm** [LR §1,§C5]. This fixes the ΔG offset / barrier-bump parameters of the Challenge-2 energy.
- **From constant-torque scans:** ramp τ and record P(open)(τ) = ⟨M/M_max⟩; the midpoint is τ_c, the steepness encodes cooperativity. The fixed-Lk τ(Lk) plateau must equal the same τ_c (ensemble consistency check).

### 3.4 Nucleation barrier and rates (Phase P4)

- **2-D metadynamics / OPES** in (Φ, M) reconstructs F(Φ,M); the saddle gives the nucleation barrier. Target the Sicard-Manghi numbers at 10-bp scale: **ΔF_op ≈ 22 kBT (opening), ΔF_cl ≈ 13–14 kBT (closure), ΔF_0 ≈ 8 kBT (formation)** [LR §7]. The barrier is of **torsional/elastic origin** — collective untwisting is rate-limiting — and **scales affinely with β·κ_φ (the melted-basin torsional modulus)**, i.e. it is a directly tunable knob via `M_melt` of the Challenge-2 energy.
- **Kramers / MFPT rate** from the 1-D PMF + position-dependent diffusion D(Φ):
  `k^{-1} = ∫_A^B dΦ e^{βF(Φ)}/D(Φ) ∫ dΦ' e^{−βF(Φ')}`  (overdamped MFPT) [LR §7; Kramers 1940],
  or the Arrhenius form `k = (ω_min ω_‡ / 2πγ) e^{−βΔF‡}`. Get D(Φ) from the autocorrelation of Φ in a restrained window, or from the string method's D(α) [LR §7; Maragliano 2006].
- **FFS cross-check** of the absolute rate using M as the order parameter (interfaces λ_0<λ_1<...): `k_AB = Φ_{A,0}·Π_i P(λ_{i+1}|λ_i)` [LR §7; Allen 2005]. FFS needs no good RC and uses natural dynamics — the gold standard for the rate, and it validates whether Φ alone suffices (compare to committor).
- **Committor / TPS spot-check** [LR §7; Dellago 1998]: harvest a handful of reactive trajectories, compute p_B at the putative transition state; if p_B≈0.5 on the Φ=Φ_‡ isosurface, Φ is a good RC; otherwise add M (2-D string).

### 3.5 Hysteresis / nonequilibrium

The finite nucleation barrier produces experimentally-observed hysteresis after torque/temperature jumps [LR §1]. Two analyses:
- **Torque/Lk ramps at varying rate** (P5): ramp τ (or turns) up and down; the loop area and the rate-dependent transition torque give the barrier via dynamic-force-spectroscopy / Bell analysis (`τ_transition` vs ln(ramp rate) → barrier and x_‡) [LR §7; Bell 1978; Dittmore 2018]. Extrapolate to zero rate for the equilibrium τ_c.
- **Two-state hopping / dwell times** at fixed τ near τ_c (analogous to AuRBT G(angle)=−kBT ln P(angle) and dwell-time kinetics [LR §1; Lebel 2014]): record residence times in each basin → rates directly, no biasing. Best done once a window straddles τ_c.

---

## 4. Concrete simulation protocol

### 4.1 System

- **Size:** N = 200 beads (≈ 2000 bp ≈ 680 nm), enough to host a 10-bp bubble (1 bead) plus elastic relaxation on both sides and to suppress finite-size barrier renormalization [LR §7; Sicard 2020]. For rate work, a shorter **N = 40–60** bead segment with one designated "weak" bead reduces cost; for sequence-dependent nucleation use a realistic genomic insert.
- **Sequence:** define per-bead types via the sequence-dependent database (`seqs_set=true` in parse_rbp.h). For controlled tests, an AT-rich bead flanked by GC-rich beads to localize nucleation [LR §3,§4; Choi 2004; Matek 2015].
- **Ends:** bead 0 clamped (reference triad); bead N−1 = handle (tension + torque/turns). Tension F = 1–2.5 pN along ẑ to stay twist-dominated [LR §1,§4].
- **Force field:** Challenge-2 double-well bond/angle/dihedral styles (log-sum-exp of two shifted Gaussians, B1) with the FENE non-extensibility add-on retained; `fix rbp/lrf` extended for unwrapped twist (§2.1).

### 4.2 Equilibration

1. **Minimize** (`min_pre_force` path already supported by the fix) to relax overlaps.
2. **Thermalize** duplex basin: NVT/Langevin, free or lightly clamped ends, 10–100 ns model time, until A and C converge (P0).
3. **Apply tension** F along ẑ, re-equilibrate; confirm Wr≈0 (§2.3).
4. **Set ensemble:** constant-τ (`fix addtorque/atom` or `fix rbp/torque`) or fixed-Lk (stiff Φ restraint). Equilibrate at the starting Φ_0 (duplex) before production.

### 4.3 Biasing schedule

- **P1 (linear response):** scan τ ∈ [0, −8] pN·nm (below τ_c); for each, ⟨Φ⟩ → slope = C_eff·2π/L. ~10 ns/point.
- **P2 (US PMF):** ~20–40 windows in Φ across the transition; per window 50–200 ns after a 10 ns window-equilibration; reweight with MBAR. Iterate window placement once using the first PMF's error profile.
- **P4 (metadynamics):** well-tempered/OPES in (Φ, M); bias factor γ≈10–20; initial hill height ~1 kBT, width ~ (0.1 turn, 0.5 step); deposit every 1–2 ps; run until the bias flattens (convergence by reweighted-FE stability and by recrossing counts). Reweight to F(Φ,M).
- **P4 (FFS):** order parameter M; interfaces at M = {0.5, 1, 2, …, M_melt}; standard direct-FFS flux + conditional probabilities.
- **P5 (hysteresis):** τ ramps at ≥4 rates spanning 2 decades; up and down; ≥10 repeats each for statistics.

### 4.4 Reweighting and analysis

- **MBAR** for all umbrella/metad data → G(Φ), F(Φ,M), with bootstrap error bars.
- **Maxwell construction** on G(Φ) → τ_c (§3.3).
- **Kramers/MFPT + D(Φ)** → rates; **FFS** → independent absolute rate; cross-check.
- **Lk = Tw + Wr** logged every frame as a correctness invariant (§2.3).

### 4.5 Observables and experimental validation targets

| Observable | How computed | Experimental target | Source |
|---|---|---|---|
| Duplex bend persistence A | tangent-tangent correlation, unbiased | 45–50 nm | [LR §1,§C3] |
| Duplex twist modulus C / C_eff | Var(Φ) or τ(σ) slope | 95–110 nm (Marko/Bryant); flag 22 nm (Gao) | [LR §1,§C3; caveat] |
| **Torque vs turns** | ⟨Φ⟩(τ) or τ(Φ*) | linear branch + plateau ("hat curve") | [LR §1,§5] |
| **Melting torque τ_c** | Maxwell on G(Φ); plateau | **−10 to −11 pN·nm** | [LR §1,§C5] |
| **P(open) vs torque** | ⟨M/M_max⟩(τ) | sigmoidal step at τ_c | [LR §1,§4] |
| Melted-basin twist | Φ in melted branch | ≈ −13 bp/turn (L-DNA) | [LR §1,§C4] |
| Nucleation barrier | metad saddle / FFS | ΔF_op≈22, ΔF_cl≈13 kBT | [LR §7] |
| Onset supercoiling | σ at first M>0.5 | σ ≈ −0.015 (F>0.5 pN) | [LR §1] |
| AT-vs-GC nucleation site | per-bead m_i | AT-rich opens first | [LR §3,§4] |
| Hysteresis loop | τ-ramp up/down | finite, rate-dependent | [LR §1] |

> **Do not** target the oxDNA ~3 pN·nm plateau (Matek 2015) — it is a model/condition-specific value, not the −10 pN·nm thermodynamic melting torque [LR §4 caveat]. Calibrate the Challenge-2 ΔG/barrier to the experimental −10 to −11 pN·nm.

---

## 5. Code-integration map (where each piece lives)

| Need | File / object | Action |
|---|---|---|
| Unwrapped per-step twist + winding | `fix_rbp_lrf.cpp/.h` | add `wind`, `tw_prev`, unwrap in `compute_lrf` Phase 2; extend grow/copy/set/pack_exchange/forward_comm; **add write_restart/read_restart** |
| Global twist CV Φ | new `compute rbp/twist` | sum unwrapped step twist over owned left atoms + Allreduce |
| Melted-step count M | new `compute rbp/melted` | sum softmax weight m_i (from double-well style) or sigmoid CV |
| Writhe / Lk check | new `compute rbp/writhe` | Klenin-Langowski O(N²) on backbone x; Fuller fast path |
| Applied torque (energy form) | new `fix rbp/torque` | add `−τ ∂Φ/∂orientation` to `atom->torque` using LRF `Jinvtp` |
| Applied torque (native) | `fix addtorque/atom` | `0 0 v_tau` on handle |
| Fixed-Lk restraint | `fix spring`/`fix restraint` on `compute rbp/twist` | `½k_Lk(Φ−Φ*)²`; read τ = ⟨k_Lk(Φ*−Φ)⟩ |
| Double-well energy + m_i export | Challenge-2 bond/angle/dihedral styles | expose softmax melted weight per junction |
| Metadynamics/OPES/reweighting | `fix plumed` | feed native CVs via per-atom property; or bias M directly; **avoid TORSION** |
| Integrator | `fix nve/dotc/langevin` (CG-DNA) or `fix nve/asphere`+`fix langevin` | rigid-body SE(3) sampling |

The existing FENE bond style (`bond_rbp_fene.cpp`) is the architectural template for adding the applied-torque scalar potential of one step coordinate with an analytic gradient: it already adds a non-quadratic scalar `E(r)` on top of the harmonic wrench and maps a scalar derivative to lab-frame force. The torque-vs-twist tilt is the rotational analogue: a scalar `E(Φ) = −τΦ` mapped to body torque via the `Jinvtp`→triad chain that `bond_rbp_fene.cpp` already executes for the wrench.

---

## 6. Pitfalls and mitigations

1. **Restart loses winding.** `fix rbp/lrf` has no `write_restart`/`read_restart`. The absolute winding integer `wind` and `tw_prev` MUST be persisted, or Φ jumps by multiples of 2π on restart and ruins reweighting. **Add restart I/O** as part of the Challenge-1 implementation; FENE already shows the restart-I/O pattern in this codebase.
2. **Per-step twist > π in one MD step.** Breaks unwrapping. Choose Δt small enough that no junction untwists by > π per step; the fix runs every step so the budget is per-MD-step, not per-dump.
3. **PLUMED TORSION is wrapped.** Never bias it for multi-turn twist. Use native unwrapped CV or the periodicity-free M [LR §7 caveat].
4. **Fuller writhe fails at antipodal points.** Switch to Klenin-Langowski near t·ẑ→−1 (plectonemes) [LR §5; Neukirch-Starostin].
5. **Ensemble mismatch.** Constant-τ and fixed-Lk must agree on τ_c via Legendre transform; if not, the bias is leaking (check the energy-form torque vs native `addtorque`) [LR §2; Moroz-Nelson].
6. **Nonequilibrium ramps overestimate τ_c.** Extrapolate τ_transition(ramp rate)→0, or measure dwell-time hopping at fixed τ near τ_c [LR §1; Dittmore 2018].
7. **Finite-size barrier renormalization.** Topologically constrained / short domains are less cooperative [LR §3,§4; Sicard 2020]; size-converge the barrier.
8. **Open-chain Lk ambiguity.** Lk strictly defined only for closed/clamped topology; for open chains use Fuller relative writhe + explicit end framing, or run minicircles for the hard invariant check [LR §5].
9. **C-modulus literature spread.** C ≈ 95–110 nm (Marko/Bryant) vs ≈22 nm (Gao-Wang) reflect force regime / bend-renormalization; pick the value matching your ensemble and tension before claiming agreement [LR §1 caveat].

---

## 7. References to literature-review sections used

§1 (experimental phenomenology, τ_melt≈−10 pN·nm, hysteresis, L-DNA basin), §2 (Maxwell/coexistence, applied-torque↔fixed-Lk Legendre, Marko 2007, Marko-Neukirch 2013), §3 (sequence-dependent ΔG, AT-rich nucleation, junction/cooperativity), §4 (CG melting models, oxDNA fixed-Lk practice, Matek 2015 caveat), §5 (Lk=Tw+Wr, Fuller/Klenin-Langowski writhe, unwrapped tracking), §6 (double-well energy forms B1/B2/B6 used for calibration), §7 (umbrella/MBAR, metadynamics/OPES, FFS, TPS/committor, string, MSM, Kramers, Bell, PLUMED caveat; Sicard-Manghi 2015 barriers; Liebl-Zacharias umbrella-under-torque; quaternion-Langevin integrators), §A (unwrapped-twist methods, dump-interval caveat), §C (parameter targets: A, C, ΔG(seq), L-DNA, τ_c).
