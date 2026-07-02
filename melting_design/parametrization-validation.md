# Parametrization & Validation Strategy: Torque-Induced DNA Melting in the CG-RBP Model

**Scope.** Concrete, implementation-ready plan to parametrize and validate a *bimodal*
(duplex / melted) sequence-dependent free energy for the existing SE(3) coarse-grained
rigid-base-pair (CG-RBP) LAMMPS force field at **~10 bp / bead** resolution, with
continuous (unwrapped) twist tracking, run in applied-torque / fixed-linking-number
ensembles, and studied as a Kramers barrier-crossing problem.

This document assumes the companion **literature review**
(`melting_design/literature_review.md`) and the design choices it recommends:
energy form **B1** (log-sum-exp of two shifted Gaussian wells), twist tracking **A1+A3**
(Y-convention per-step deformation + winding counter, cross-checked by Lk = Tw + Wr).

All energies are reported in **kBT at T = 300 K** unless noted.
Conversions used throughout: 1 kcal/mol = 1.688 kBT @ 310 K = 1.677 kBT @ 300 K;
kBT(300 K) = 4.114 pN·nm. A "bead" spans **n_bp = 10** base-pair steps (one helical turn).

---

## 0. Code anchors (what we reuse, what we add)

Grounding the design in the actual source under `src/CG-RBP/`:

| Concern | Existing code fact | What we do |
|---|---|---|
| Per-step deformation `Y = (Ω_d, w_d)` | `bond_rbp_fene.cpp` computes `Φ_d = log(Sᵀ R)` (Y conv) or `Ω - srot` (X conv); energy `½ YᵀMY`; `Yv[6]` assembled at line ~396 | Add a **second basin** `(Y0_melt, M_melt)` and blend via log-sum-exp; reuse the *same* gradient-to-torque path (`leftJacobianInverseTransposed`, `Smat`) |
| Twist component | `so3::rotmat2euler` caps `|Ω| = acos((trR−1)/2) ∈ [0,π]` (so3.h:98–155); branch-cut handling near θ=π | Keep for in-frame `Φ_d`; add **winding counter** `nw` so reported twist = `Φ_d,3 + 2π·nw` |
| Persistent per-junction state | `fix_rbp_lrf` holds `triads/fwd_euler/fwd_Jinvtp` keyed to **left atom**, with `grow/copy/set_arrays`, `pack_exchange`, `pack_forward_comm`; **no** `write_restart/read_restart` | Add `int *nwind` per left-atom; extend the same callbacks; **add** restart I/O (winding must survive restarts) |
| Non-quadratic scalar add-on | `bond_rbp_fene.cpp` adds one-sided FENE `E(r)` on top of wrench, analytic `dE/dr` along bond | Architectural template for the **barrier bump** `B·exp(−(Tw−Tw*)²/2w²)` on the twist coordinate |
| Coefficient format | `parse_rbp.h`: `Ystatic[6]` + 6×6 `M` (48 coeffs/type); flags `subtract_groundstate, seqs_set, closed, unit_length` | Add per-type **melted block** (`Ystatic_melt[6]`, `M_melt` 21, `ΔG`, `B`, `Tw*`, `w`, `ρ`); add metadata flag `bimodal` |

The existing model already supplies the **duplex basin** `(Y0_dup(seq), M_dup(seq))`
per bond/angle/dihedral type. The new work is the melted basin, the offset, the barrier,
and the winding bookkeeping.

---

## 1. Duplex basin (stiff well) — reuse + verify

### 1.1 Source of stiffness and ground state

Reuse the existing CG-RBP coarse-grained parameters. Their provenance and the
recommended reference chain:

- **Microscopic stiffness:** cgDNA / cgDNA+ / cgNA+ (Gonzalez–Petkevičiūtė–Maddocks 2013;
  Sharma–Patelli–Maddocks 2023) provide sequence-dependent rigid-base ground states and
  6×6 (rigid-base) stiffness from extensive MD; Olson 1998 (knowledge-based) and
  Lankaš 2003 (MD covariance, `K = kBT·cov⁻¹`) are independent cross-checks for the
  per-dinucleotide intrinsic twist (~34–36°/step) and force constants.
- **Coarse-graining to 10 bp:** PolyCG (Skoruppa–Schiessel 2025) — the *same* framework
  underlying this model — marginalizes the per-step Gaussian to a per-bead Gaussian by
  composing junctions and propagating the 6×6 covariance. This **automatically** sets the
  10-bp intrinsic twist (≈ +1 turn/bead) and the coarse stiffness.

### 1.2 How 10-bp coarse-graining fixes intrinsic twist & stiffness

Intrinsic twist per bead (the third component of `Ystatic`/`srot`):

```
Ω₃⁰_bead(seq) = Σ_{k=1..n_bp} twist_step_k(seq)   (mod nothing — accumulate, see §4)
              ≈ 10 × 34.3°  ≈ 343° ≈ 5.99 rad ≈ +1.0 turn
```

With ψ_B = 0.598 rad/bp (Marko–Neukirch 2013), one 10-bp bead is **+5.98 rad ≈ +0.95
turn** of intrinsic twist. *Note:* this exceeds the π cap of `rotmat2euler`; the duplex
ground state is stored in `Ystatic` directly (never round-tripped through a rotation
matrix), and the Y-convention factors `S` out so the *deformation* `Φ_d` stays small.
This is precisely why the existing Y convention is mandatory and why winding tracking
(§4) is needed once melting shifts the live twist by ~−1 turn.

Coarse stiffness (length-scale-dependent apparent moduli; Becker–Everaers 2007):

| Quantity | Per-step (cgDNA-level) | 10-bp bead (target after CG) | Source |
|---|---|---|---|
| Bend persistence A | — | **45–50 nm** | Marko–Siggia; Skoruppa 2017 |
| Twist modulus C | — | **95–110 nm** (intrinsic) / ~22–80 nm (C_eff, renormalized) | Marko–Neukirch (95); Bryant (~105); Gao–Wang (22) |
| Twist–bend G | — | **30 nm** (grooved) | Skoruppa 2017; Nomidis 2019 |

`M_dup` (6×6) for a bead is obtained from A, C, G via the standard TWLC map
(Nomidis–Skoruppa–Carlon 2019), divided by bead length a = n_bp × 0.34 nm = 3.4 nm:

```
βE_dup,step = (1/2a)[ A₁ ω₁² + A₂ ω₂² + C ω₃² + 2G ω₂ω₃ ]   (rotational block)
```

**Action item D1.** Verify the *existing* coarse `M_dup` reproduces A≈45 nm, C≈100 nm,
G≈30 nm by the fluctuation route `M⁻¹ = ⟨ΔY ΔYᵀ⟩` (Skoruppa 2017 covariance fit) in a
free (zero-torque) simulation of a 50-bead duplex. This is a *consistency check on the
unchanged model*, and the regression baseline for all later work.

---

## 2. Melted basin (soft well)

### 2.1 Target parameters (L-DNA, the experimentally defined melted phase)

The melted basin is **not** a structureless open bubble: under torque it is the
left-handed L-DNA phase with its own intrinsic twist and (much softer) elasticity
(Sheinin–Wang 2011; Marko–Neukirch 2013). Use L-DNA as the primary target, with a
fully-denatured ssDNA limit as a sensitivity bracket.

| Quantity | Duplex (B) | **Melted (L-DNA), primary** | ssDNA limit (bracket) | Source |
|---|---|---|---|---|
| Intrinsic twist | +34.3°/bp (+0.598 rad/bp) | **−13 bp/turn = −0.393 rad/bp** | ≈ 0 | Sheinin 2011; Marko–Neukirch 2013 |
| Twist shift / 10-bp bead | (ref) | **ΔTw ≈ −2π (one full turn)** | ≈ −1 turn | derived below |
| Bend persistence A | 45–50 nm | **3–7 nm** | 0.7–1.5 nm | Sheinin (3); Marko–Neukirch (7); Maffeo (ss) |
| Twist modulus C | 95–110 nm | **19–20 nm** | ≈ 0 | Sheinin; Marko–Neukirch |
| Rise / bp | 0.34 nm | **0.459–0.48 nm** | ~0.5–0.7 nm | Marko–Neukirch; Smith 1996 |

### 2.2 Twist shift of the melted ground state (the −1 turn)

Per 10-bp bead, the difference in integrated intrinsic twist between basins:

```
ΔΩ₃ = n_bp × (ψ_L − ψ_B)
    = 10 × (−0.393 − 0.598) rad
    = 10 × (−0.991) = −9.91 rad ≈ −1.58 turn
```

Per *contact*, B-DNA stores +0.95 turn and L-DNA −0.62 turn, so a 10-bp segment that
melts releases **≈ −1.58 turns of twist** into the rest of the molecule (writhe + bulk
twist). Matek 2015 (oxDNA) reports a melted bubble dumps ~1 turn of undertwist locally,
consistent at this resolution. **Set the melted ground state**:

```
Ystatic_melt[2] (twist)  = Ystatic_dup[2] + n_bp·(ψ_L − ψ_B)   ≈ Ystatic_dup[2] − 9.91 rad
Ystatic_melt[0,1]        = Ystatic_dup[0,1]                      (tilt/roll intrinsic ≈ 0)
Ystatic_melt[3,4,5] (w)  : rise scaled 0.34→0.459 nm  ⇒  |w_melt| = (0.459/0.34)|w_dup| ≈ 1.35|w_dup|
```

**Storage subtlety.** `Ystatic_melt[2]` exceeds π and must **not** be round-tripped
through `rotmat2euler`. Define the melted `Smat_melt = euler2rotmat(srot_melt)` once at
coefficient-read time (the Rodrigues form `euler2rotmat` in so3.h:50 is valid for any
magnitude). The deformation `Φ_d,melt = log(Smat_meltᵀ R)` then stays small near the
melted minimum — exactly mirroring the duplex Y convention. The −1-turn *jump between
basins* is carried by the winding counter (§4), not by either basin's `Φ_d`.

### 2.3 Melted stiffness matrix

```
A_L = 3–7 nm,  C_L = 19–20 nm,  G_L ≈ 0  (no grooves in the melted state)
M_melt rotational block (per bead, a_L = n_bp·0.459 nm):
   M_melt = diag(A_L, A_L, C_L)/a_L     (then convert to kBT units)
Translational block: soften ~3–5× vs duplex (ssDNA stretch compliance);
   keep FENE non-extensibility (bond_rbp_fene) active so the soft basin cannot blow up.
```

The **ratio** C_dup/C_melt ≈ 100/20 = 5 (L-DNA) up to ≈ ∞ (ssDNA) is the single most
important number: it is the origin of the cooperativity and a primary barrier knob
(Sicard–Manghi 2015: barrier ∝ β·κ_φ of the duplex). Start with C_melt = 20 nm and treat
it as a tunable in §5.

---

## 3. Sequence-dependent inter-basin offset ΔG(seq)

### 3.1 Per-bead duplex stabilization from nearest-neighbor parameters

`ΔG(seq)` is the free energy by which the **duplex basin sits below the melted basin** at
zero torque. Build it from SantaLucia unified nearest-neighbor (NN) parameters
(SantaLucia 1998; SantaLucia–Hicks 2004), summed over the n_bp steps in each bead:

```
ΔG_bead(seq) = ΔG_init + Σ_{steps∈bead} ΔG_NN(step)   (+ termAT for terminal beads)
```

ΔG37 per step (kcal/mol, 1 M Na⁺):

| step | ΔG37 | step | ΔG37 |
|---|---|---|---|
| GC/GC | −2.24 | GA/CT | −1.30 |
| CG/CG | −2.17 | CT/GA | −1.28 |
| GG/CC | −1.84 | AA/TT | −1.00 |
| CA/GT | −1.45 | AT/AT | −0.88 |
| GT/CA | −1.44 | TA/TA | −0.58 |

Initiation +0.98 (G·C end) / +1.03 (A·T end); symmetry +0.43.

**Per-bead range** (9 internal steps): AT-rich ≈ 9×(−0.7) ≈ −6.3 kcal/mol ≈ **−10.6 kBT**;
GC-rich ≈ 9×(−2.2) ≈ −19.8 kcal/mol ≈ **−33 kBT**. So at 300 K, `ΔG(seq)` ranges
**≈ 7–24 kBT** per bead favoring duplex (consistent with Benham per-bp 1.0–3.5 kcal/mol;
Cocco–Marko ΔG_AT ≈ 1.1 kBT/bp, ΔG_GC ≈ 3.5 kBT/bp).

### 3.2 Temperature and salt corrections (do NOT skip)

Torque tweezers run at ~100–150 mM Na⁺, not 1 M. Apply:

```
ΔG(T) = ΔH − T·ΔS                         (full ΔH/ΔS tables, SantaLucia 1998)
ΔG_salt = ΔG(1M) + 0.114·N_bp·ln([Na⁺])   per duplex (Owczarzy 2004, monovalent)
                                           (Mg²⁺/mixed: Owczarzy 2008, R=√[Mg]/[Mon])
```

AT-rich beads are more salt-sensitive — the correction is per-bp, so a 10-bp bead at
150 mM is destabilized by ≈ 10×0.114×ln(0.15) ≈ −2.2 kcal/mol ≈ **−3.7 kBT** relative to
1 M, i.e. ΔG shrinks. Cross-check the bead-averaged value against Marko–Neukirch
ε_M = 2.7 kBT + 0.2·ln([Na⁺]/150mM) per bp ⇒ ~27 kBT/bead averaged — note this is the
*formation* energy of L-DNA which already nets the L-DNA elastic ground-state energy; see
§3.4 for the bookkeeping that avoids double counting.

### 3.3 Mapping per-dinucleotide ΔG to a per-bead offset

Each bead = one bond/angle/dihedral *junction* keyed to its left atom (the existing
`fix_rbp_lrf` convention). Because the database is already per-type and per-junction
(`parse_rbp.h`), the natural mapping is:

```
for each junction type t (covering bp index i..i+n_bp):
    ΔG_t = ΔG_init·[t is terminal] + Σ_{k=i..i+n_bp-1} ΔG_NN(seq[k:k+2])
            + salt_correction(n_bp, [Na+]) + (ΔH−TΔS adjustment to sim T)
    store ΔG_t in the melted block of the coefficient record
```

This is a one-time **preprocessing step in `cgRBPTools`** (the Python side): read the
sequence, slide a 10-bp window, emit `ΔG_t` per type into the database file alongside
`Ystatic_melt` and `M_melt`. The C++ side only reads it.

### 3.4 Bookkeeping: separate elastic energy from the offset (avoid double counting)

The log-sum-exp blend (B1, §5) already contains the *elastic* energy of each basin. The
offset `ΔG(seq)` must therefore be the **pure basin-bottom free-energy difference**:

```
E_blend = −(1/ρ) ln[ exp(−ρ E_dup) + exp(−ρ (E_melt + ΔG_offset)) ]
ΔG_offset(seq) = ΔG_NN(seq)            (basin-bottom duplex stabilization, §3.1–3.3)
```

Crucially `E_dup`, `E_melt` are each ≈ 0 at their own minima, so `ΔG_offset` is exactly
the well-bottom gap. The L-DNA *elastic ground-state* energy relative to B (the cost of
sitting in a softer, retwisted well) is captured by `M_melt` and `Ystatic_melt`, **not**
by `ΔG_offset`. When calibrating against Marko's ε_M (a combined number), back out the
elastic part analytically before assigning `ΔG_offset` (§5.3 calibration).

---

## 4. Continuous (unwrapped) twist tracking

This is **Hard Problem (1)** and a prerequisite for everything in §5–§7. Recommended:
**A1 (per-step deformation + winding counter)** as the working twist; **A3 (Lk=Tw+Wr)**
as a per-frame topology cross-check / collective variable.

### 4.1 Per-junction winding counter (A1)

Attach an integer `nwind` to each junction (keyed to left atom, exactly like
`fwd_euler`). At each `pre_force`, after computing the bare per-step deformation twist
`φ = Φ_d[2]` (already small in the Y convention), unwrap against the previous frame:

```cpp
// in FixRBPLRF: new per-atom int array nwind[i], double prev_phi[i]
double phi = fwd_euler[i][2];                 // current bare twist deformation (Y conv)
double dphi = phi - prev_phi[i];
if      (dphi >  M_PI) nwind[i] -= 1;          // wrapped +π→−π : true twist decreased by 2π
else if (dphi < -M_PI) nwind[i] += 1;
prev_phi[i] = phi;
// unwrapped per-junction twist available to energy/CV:
double twist_unwrapped = phi + 2.0*M_PI*nwind[i];
```

**Total chain twist** (collective variable for §6–§7):
`Tw = Σ_junctions (Φ_d[2] + 2π·nwind)`.

Validity condition: the per-junction twist must change by < π between consecutive
`pre_force` calls. At 10-bp resolution and a sane timestep this holds easily; enforce it
with an assertion that flags `|dphi| > 0.8π` (should never fire in equilibrium dynamics).

### 4.2 Required fix_rbp_lrf changes (concrete)

`fix_rbp_lrf` currently has **no restart I/O** — winding number is path-dependent state
that *must* persist. Additions:

```cpp
// new persistent per-atom state
int    *nwind;        // winding counter per left-atom junction
double *prev_phi;     // last-frame bare twist deformation, for unwrap

// extend existing callbacks (mirror triads/fwd_euler handling):
grow_arrays:     memory->grow(nwind,...); memory->grow(prev_phi,...);
copy_arrays:     copy nwind[i]->nwind[j], prev_phi
set_arrays:      nwind[i]=0; prev_phi[i]=0
pack_exchange / unpack_exchange:  +2 doubles (cast int)
pack_forward_comm:                ghosts need unwrapped twist for coupled (3-/4-pt) terms

// NEW (does not currently exist — add it):
void write_restart(FILE*);  // dump nwind, prev_phi per owned atom
void read_restart(FILE*);   // restore; without this, restarts reset winding to 0
int  size_restart(int);  void *extract(...);  // per-atom restart hooks
```

Mirror the restart pattern already shown in `bond_rbp_fene.cpp::write_restart/read_restart`
(per-type) but at the **per-atom** level (use `atom->add_callback(Atom::RESTART)` and the
`pack_restart/unpack_restart` Fix hooks).

### 4.3 Topology cross-check (A3)

Each dump frame, compute writhe Wr from bead centerlines (Klenin–Langowski exact
solid-angle sum; Fuller single-integral vs straight reference for speed in extended
configs) and verify `Lk = Tw + Wr` is conserved (integer, fixed by the ensemble). This is
an *analysis-time* check (Python in `cgRBPTools`), not in the force loop. It catches
unwrap failures and is the natural CV for the fixed-Lk ensemble. Fall back to the full
double sum when `t·t₀ → −1` (antipodal/plectoneme; Neukirch–Starostin 2008).

### 4.4 PLUMED caveat

PLUMED's `TORSION` CV is wrapped to (−π,π]. Any biasing of twist (§7) must feed the
**unwrapped** `Tw` via a `fix colvars`/custom variable, **or** bias the periodicity-free
**melted-step count** `m = Σ s_i` (§6, B3) instead.

---

## 5. Bimodal energy & barrier (the new force law)

### 5.1 Primary form B1 (log-sum-exp of two shifted Gaussians) + barrier bump

Per junction (bond/angle/dihedral), with `Y = (Φ_d, w_d)`:

```
E_dup(Y)  = ½ (Y − 0)ᵀ M_dup  (Y)              [Y already deformation from S_dup; min at 0]
E_melt(Y) = ½ (Y − Y0Δ)ᵀ M_melt (Y − Y0Δ)      [Y0Δ = melted min in the dup-deformation frame]
E(Y)      = −(1/ρ) ln[ exp(−ρ E_dup) + exp(−ρ (E_melt + ΔG_offset)) ]
                 + B·exp( −(Tw_loc − Tw*)² / (2 w²) )           [explicit barrier bump]
```

Force/torque (analytic, cheap):

```
p_dup  = exp(−ρ E_dup) / Z,   p_melt = 1 − p_dup,   Z = exp(−ρE_dup)+exp(−ρ(E_melt+ΔG))
∂E/∂Y  = p_dup·∂E_dup/∂Y + p_melt·∂E_melt/∂Y      (softmax-weighted basin gradients)
```

Each `∂E_{basin}/∂Y` is computed **exactly as the existing single-well code already does**
(the `A = Mr·Ω_d + Mtr·w_d`, `B = ...` blocks in `bond_rbp_fene.cpp`), once per basin,
then weighted by `p_dup`/`p_melt`. The torque mapping through `leftJacobianInverseTransposed`
and `Smat` is reused verbatim per basin. The barrier-bump term is a scalar function of the
unwrapped local twist `Tw_loc`; its gradient w.r.t. the twist DOF is added along the twist
torque axis exactly like FENE adds `dE/dr` along the bond (the `bond_rbp_fene` template).

### 5.2 Reduced 1D form B2 (for Kramers / PMF validation only)

```
V(φ) = A₄φ⁴ + A₂φ² + A₁φ,   φ = Tw − Tw0,  A₂<0
minima ±√(−A₂/2A₄);  barrier = A₂²/(4A₄);  A₁ = −τ·(dTw/dφ lever)  (applied torque tilt)
```

Used to fit the metadynamics PMF and extract a Kramers rate; not the MD force law.

### 5.3 Parameter calibration (fixing the four knobs: ΔG_offset, Y0_melt, M_melt, B, ρ)

Calibrate against the **macroscopic constraint** (B6, Marko 2007 Maxwell construction):
coexistence (degenerate basins under torque) must occur at the experimental melting torque.

```
Coexistence torque:  τ_c = (E_melt_min + ΔG_offset − E_dup_min) / (Tw0_dup − Tw0_melt)
Target:              τ_c = −10 to −11 pN·nm  (Sheinin 2011; Bryant 2003; Mosconi 2009; Marko 2007)
```

With `Tw0_dup − Tw0_melt ≈ +9.91 rad/bead` (§2.2) and `τ_c·ΔTw` giving the tilt at
coexistence, solve for the basin-bottom gap that, *combined with* `ΔG_offset(seq)`,
reproduces −10 pN·nm for the sequence-averaged bead. This anchors the *average*; the
sequence dependence is then carried by `ΔG_offset(seq)` from §3, giving AT-rich beads a
lower melting torque (melt first) — matching Sheinin 2011 (strongly sequence-dependent
below 5 pN) and Liebl–Zacharias (TATA melts first).

**Do NOT** target the oxDNA ~3 pN·nm plateau (Matek 2015) — that is a model/condition
artifact, not the thermodynamic melting torque.

### 5.4 Barrier height calibration (B, ρ, and C_melt)

Three coupled knobs set the kinetic barrier; calibrate against bubble nucleation cost:

| Knob | Effect | Target |
|---|---|---|
| `ρ` (softmax sharpness) | ρ→∞ sharp cusp/high barrier; small ρ merges wells | coarse barrier shape |
| `B`, `w` (bump) | independent, localized barrier at `Tw*` (crossover twist) | fine kinetic tuning |
| `C_melt` (melted twist stiffness) | barrier ∝ β·κ_φ of stiff side (Sicard–Manghi) | physical origin |

**Targets** (10-bp bubble): nucleation barrier **ΔF_op ≈ 22 kBT**, closure
**ΔF_cl ≈ 13–14 kBT**, formation **ΔF₀ ≈ 8 kBT** (Sicard–Destainville–Manghi 2015;
Manghi–Destainville 2016). The junction/cooperativity cost (Benham a ≈ 10–11 kcal/mol
≈ 17–18 kBT; Poland–Scheraga σ_coop ~ 10⁻⁴–10⁻⁵) is reproduced *either* by the bump `B`
(B1) *or* by an explicit inter-bead junction penalty `J` (form B3) if multi-bead
cooperativity must be modeled. Buckling-only barriers are smaller (~10 kBT, Dittmore 2018)
and provide a low-end sanity bracket.

### 5.5 Optional cooperativity (B3, only if needed)

If single-junction B1 under-predicts cooperativity (melting bubbles too short / too
frequent), add a discrete spin `s_i ∈ {dup, melt}` per junction and an inter-bead
domain-wall penalty `J` between unlike neighbors (Storm–Nelson / Benham transfer matrix).
At 10-bp resolution this may be unnecessary (one bead ≈ one nucleation unit); defer until
validation ladder rung (c) demands it.

### 5.6 Applied-torque / fixed-Lk ensembles

- **Constant torque:** add `−τ·Tw` (Gibbs/Legendre) via `fix addtorque`/`fix addtorque/atom`
  on terminal beads (clamp the other end with `fix rigid` or angular springs). Direct,
  matches optical/magnetic torque-wrench experiments.
- **Fixed Lk:** clamp both ends' orientation (impose target integer turns), measure the
  conjugate torque from `⟨Tw⟩` fluctuations or reaction torque on the clamp. Matches
  magnetic-tweezers `ΔLk = n` (Strick 1996). Use `Lk = Tw + Wr` (§4.3) as the constraint
  CV. Integrator: quaternion-Langevin `fix nve/dotc/langevin` (Davidchack 2017) /
  `fix nve/asphere` — note these advance SO(3) and do **not** store winding, so §4.1 must
  reconstruct it.

---

## 6. Collective variables

| CV | Definition | Periodicity-safe? | Use |
|---|---|---|---|
| `Tw` (unwrapped) | `Σ (Φ_d[2] + 2π·nwind)` | yes (via §4) | torque ensembles, Kramers RC, metad |
| `Lk` | `Tw + Wr` (Klenin–Langowski) | yes | fixed-Lk constraint, topology check |
| `m` (melted-step count) | `Σ s_i`, `s_i = sigmoid(ρ(E_dup−E_melt−ΔG))_i` | yes (discrete) | FFS order parameter, bubble size |
| `Wr` | Gauss double integral / Fuller | yes | plectoneme vs bubble partition |

Twist alone is the *natural* RC (Sicard–Manghi: collective untwisting is rate-limiting;
Dasanna 2013: winding is the slow coordinate). Use `m` for FFS (periodicity-free) and
validate twist-as-RC via committor / transition-path sampling (Dellago 1998).

---

## 7. Validation ladder

A staged ladder, each rung with concrete target numbers and a comparison dataset. **Do not
proceed to the next rung until the current one passes.** Regression each rung in CI.

### Rung (a) — Single junction: two-state occupancy vs torque  *[unit test]*

**System:** one bond (two beads), one junction type, fixed sequence, applied constant
torque scan `τ ∈ [−20, +5] pN·nm`. Long Langevin run per τ; histogram `Tw_loc`.

**Observables / targets:**
- Bimodal `P(Tw_loc)` with two peaks at `Tw0_dup` and `Tw0_dup − 9.91 rad`.
- Occupancy `p_melt(τ)` is a sigmoid crossing 0.5 at `τ_c`.
- **`τ_c = −10 to −11 pN·nm`** (sequence-averaged bead) — Sheinin 2011 / Marko 2007.
- Sequence shift: AT-rich `τ_c` higher (less negative) than GC-rich by the amount implied
  by `ΔG_offset` difference / ΔTw.
- PMF `−kBT ln P(Tw_loc)` fits B2 quartic; barrier matches the input `B`/`ρ`/`C_melt`
  setting to within the Kramers prefactor.

**Pass criterion:** `τ_c` within ±1 pN·nm of −10; basin minima within 5% of design;
detailed-balance check (forward/back rates from dwell times give the same ΔG as
−kBT ln(p_melt/p_dup)).

### Rung (b) — Short duplex: torque–turns & extension–turns  *[10–50 beads]*

**System:** 100–500 bp (10–50 beads), one end clamped, other end under `fix addtorque`
(constant-torque) **and** a "hat curve" run at fixed Lk (count turns, Strick 1996 geometry).
Stretching tension F ∈ {0.3, 0.5, 1, 3} pN.

**Observables / targets:**
- **Torque vs ΔLk:** linear elastic branch (slope set by C_eff), then a **flat plateau at
  τ_c ≈ −10 pN·nm** on underwinding (constant-torque coexistence = the double-well
  fingerprint; Marko 2007; Mosconi 2009).
- **Extension vs turns ("hat curve"):** symmetric at low F; on underwinding above
  F_c ≈ 0.5 pN the negative-turn arm **plateaus** (denaturation absorbs turns instead of
  writhe; Strick 1998).
- Onset of melting near **σ ≈ −0.015** (above 0.5 pN); melted fraction grows linearly with
  further underwinding (lever rule).
- Elastic branch slope reproduces C_eff (Moroz–Nelson renormalization; pick the value
  matching the ensemble/force — flagged in lit review).

**Compare to:** Strick 1996/1998 hat curves; Mosconi 2009 / Sheinin 2011 torque plateaus;
Marko 2007 analytic Maxwell construction.

**Pass criterion:** plateau torque −10±1.5 pN·nm; F_c for the underwinding plateau in
0.4–0.7 pN; backbending (negative torsional compliance) visible near coexistence
(Oberstrass 2013).

### Rung (c) — Bubble nucleation statistics & sequence dependence  *[enhanced sampling]*

**System:** designed sequences with AT-rich and GC-rich blocks (mimic Oberstrass 2012
designed substrates), fixed σ ≈ −0.03 to −0.06. Well-tempered metadynamics on
(`Tw` or `m`, optionally 2D with `Wr`) via PLUMED; or FFS on `m`.

**Observables / targets:**
- Bubbles **nucleate in AT-rich beads** (Matek 2015: 84–91% AT in tip bubbles;
  Choi 2004; SIDD profiles).
- 2D free-energy landscape `G(Tw, m)` shows two basins + saddle; **nucleation barrier
  ΔF_op ≈ 22 kBT, closure ≈ 13 kBT** for a single 10-bp bubble (Sicard–Manghi 2015).
- Per-bead melting probability profile correlates with `ΔG_offset(seq)` / SIDD G(x).
- Kramers rate `k = (ω_min ω_b / 2πγ) exp(−βΔF)` from the PMF + diffusion `D(s)`
  (string method / MFPT); opening times reach **ms**, closure **µs–ms** (Sicard 2015/2020).
- **Hysteresis** after torque jumps (finite barrier) reproduced (Strick; Argudo–Purohit).

**Compare to:** Oberstrass 2012/2013 (per-bp ΔG + junction energy, constant-torque
plateaus); Sicard–Manghi 2015/2020 barriers and rates; SIDD/WebSIDD per-bp profiles.

**Pass criterion:** correct AT-first localization; barrier within ~30% of 22 kBT (tunable
via §5.4); rate order-of-magnitude (ms opening) correct.

### Rung (d) — Cross-model & experimental benchmark  *[capstone]*

**Compare against:**
- **oxDNA / oxDNA2** (LAMMPS CG-DNA, Henrich 2018) on a *matched* short construct: Lk=Tw+Wr
  partition, bubble localization, plectoneme/bubble competition (F threshold ~2.5 pN,
  Matek 2015). Note oxDNA's ~3 pN·nm plateau is model-specific — compare *mechanism*
  (where/how bubbles nucleate), not the plateau magnitude.
- **Experiment:** torque plateau −10 to −11 pN·nm (Sheinin 2011, Bryant 2003, Mosconi 2009);
  L-DNA elastic constants (Sheinin 2011); melting onset σ ≈ −0.015 and completion σ ≈ −1.8
  (Allemand 1998); sequence-resolved transitions (Oberstrass 2012/2013); all-atom G(σ)
  three-regime shape and TATA-first melting (Liebl–Zacharias).
- **NN thermodynamics:** finite-duplex `Tm`/ΔG from the bead `ΔG_offset(seq)` vs SantaLucia
  predictions (zero-torque consistency).

**Pass criterion:** quantitative agreement on τ_c, L-DNA moduli, melting-onset σ, and
sequence ordering of melting; qualitative agreement with oxDNA mechanism and all-atom
G(σ) shape.

---

## 8. Prioritized implementation plan

| Priority | Task | Files | Depends on |
|---|---|---|---|
| **P0** | Winding counter + restart I/O in `fix_rbp_lrf` | `fix_rbp_lrf.{h,cpp}` | — |
| **P0** | Verify unchanged duplex `M_dup` ⇒ A,C,G (D1) | analysis (`cgRBPTools`) | — |
| **P1** | Extend coeff record: melted block + `ΔG,B,Tw*,w,ρ`; metadata `bimodal` | `parse_rbp.h`, `cgRBPTools` | P0 |
| **P1** | Melted `Smat_melt` from >π twist (reuse `euler2rotmat`) | `bond_rbp.cpp` etc. | P1 coeff |
| **P1** | Log-sum-exp blend in bond/angle/dihedral `compute()` (softmax-weighted dual-basin gradient) | `bond_rbp*.cpp`, `angle_rbp.cpp`, `dihedral_rbp.cpp` | P1 |
| **P2** | Barrier bump on unwrapped twist (FENE-style scalar add-on) | `bond_rbp*.cpp` | P1, P0 winding |
| **P2** | NN→ΔG(seq) preprocessor (SantaLucia + salt/T) | `cgRBPTools` | — |
| **P2** | Calibrate τ_c=−10, barrier=22 kBT (rungs a,b) | input scripts | P1,P2 |
| **P3** | Lk=Tw+Wr analysis (Klenin–Langowski) + PLUMED custom-CV | `cgRBPTools`, PLUMED | P0 |
| **P3** | Metadynamics/FFS, Kramers rates (rung c); oxDNA/expt benchmark (rung d) | input scripts | all |
| **P4** | Optional B3 inter-bead cooperativity `J` | new angle-like term | only if rung c fails |

**Critical path:** P0 winding (with restart) → P1 dual-basin energy → P2 calibration →
rungs (a),(b) → P3 rates → rung (d). The winding counter is the long pole: without it
neither the −1-turn basin shift nor the twist CV is meaningful, and it is the one piece
with no existing template (must add per-atom restart hooks `fix_rbp_lrf` lacks today).

---

## 9. Parameter table (starting values, sequence-averaged 10-bp bead, 300 K, 150 mM Na⁺)

| Symbol | Meaning | Value | Source |
|---|---|---|---|
| n_bp | bp/bead | 10 | design |
| a_dup | bead length (B) | 3.40 nm | 10×0.34 |
| a_melt | bead length (L) | 4.59 nm | 10×0.459 |
| Ω₃⁰_dup | duplex intrinsic twist | +5.98 rad (+0.95 turn) | ψ_B=0.598 |
| Ω₃⁰_melt | melted intrinsic twist | −3.93 rad (−0.62 turn) | ψ_L=−0.393 |
| ΔTw | basin twist shift | −9.91 rad (−1.58 turn) | §2.2 |
| A_dup, C_dup, G_dup | duplex moduli | 45, 100, 30 nm | Marko–Neukirch; Skoruppa 2017 |
| A_melt, C_melt, G_melt | melted moduli | 5, 20, 0 nm | Sheinin 2011 |
| ΔG_offset(seq) | basin-bottom gap | 7–24 kBT (AT→GC) | SantaLucia + salt |
| τ_c | melting torque (target) | −10 to −11 pN·nm | Sheinin; Marko 2007 |
| ΔF_op | nucleation barrier (target) | ≈ 22 kBT | Sicard–Manghi 2015 |
| ΔF_cl | closure barrier (target) | ≈ 13 kBT | Sicard–Manghi 2015 |
| σ_onset | melting onset | ≈ −0.015 | Allemand 1998; Strick 1998 |
| F_c | underwinding-plateau force | ≈ 0.5 pN | Strick 1998 |

---

## 10. Risks & flagged uncertainties

- **Torsional modulus spread:** C ≈ 95–110 nm (Marko–Neukirch/Bryant, intrinsic) vs
  ~22 nm (Gao–Wang, extended C_eff). Pick the value matching the *ensemble* (intrinsic for
  the bare `M_dup`; C_eff emerges from bend fluctuations in the simulation). Validate via D1.
- **Melting-torque magnitude:** robustly −10 to −11 pN·nm in tweezers; a higher critical
  torque (~−31 pN·nm) appears in some elastic-rod analyses (Vologodskii–Frank-Kamenetskii).
  Target −10; treat −31 as an upper sensitivity bracket.
- **oxDNA plateau (~3 pN·nm)** is model/condition-specific — use for *mechanism*
  comparison only, never as the τ_c target.
- **ΔG double-counting** (§3.4): the offset must be the basin-bottom gap, elastic
  ground-state difference goes in `M_melt`/`Ystatic_melt`. Verify by checking E_blend at
  each basin minimum equals 0 and ΔG_offset respectively.
- **Unwrap fragility** in tightly writhed configs (per-step twist > π/frame): mitigate with
  small dump interval + the `|dphi|>0.8π` assertion; A2 (parallel transport) is the fallback.
- **Restart correctness:** winding is path state — a restart without `read_restart` silently
  resets it, corrupting any long fixed-Lk run. This is the highest-severity code risk; test
  explicitly (write restart mid-melting, reload, confirm `Tw` continuous).
