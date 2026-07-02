# Tracking Unwrapped Twist / Winding Number in the CG-RBP SE(3) Model

**Scope.** Concrete implementation design for tracking *continuous, unwrapped* twist
(winding number) per step in the LAMMPS CG-RBP code, as required by torque-induced
melting at ~10 bp/bead (a ~1-turn untwist event). This document covers (A) why the
current representation breaks past ±180°, (B) a robust continuity-based unwrapping
scheme with a periodicity-free twist extractor, (C) where and how to store the state
(`fix_rbp_lrf` persistent arrays + restart + MPI/ghost), (D) exposing it as a
collective variable / compute, (E) force/torque consistency, (F) edge cases, and
(G) a step-by-step plan with validation tests.

Grounded in the actual source:
`fix_rbp_lrf.{h,cpp}`, `so3.h`, `bond_rbp.cpp`, `bond_rbp_fene.cpp`,
`dihedral_rbp.cpp`, `angle_rbp.cpp`, `parse_rbp.h`.

Literature anchors (see master review): Brackley–Morozov–Marenduzzo 2014 (LAMMPS
unwrapped twist sum Σ(α_i+γ_i)); Bergou et al. 2008 (Discrete Elastic Rods,
parallel transport + scalar twist); White 1969 / Fuller 1971,1978 / Klenin–Langowski
2000 (Lk = Tw + Wr); Dennis–Hannay 2005 (integer-jump structure of twist/writhe);
Sicard–Destainville–Manghi 2015 and Dasanna et al. 2013 (untwisting is the
rate-limiting reaction coordinate; CV = total twist Φ = Σφ_i).

---

## 0. Notation and the existing data path

Per step *i* the fix computes (`fix_rbp_lrf.cpp`, `compute_lrf`, Phase 2):

- Triads `T_i = quat_to_mat(q_i)` for all atoms incl. ghosts (Phase 1).
- For each bond `(id1, id2, btype)` keyed to the **left atom `id1`**:
  - `R = T1^T T2`  (`lamath::mul_AtB(T1,T2,R)`)
  - **X convention** (`subtract_groundstate=true`): `Omega = rotmat2euler(R)`, store in
    `fwd_euler[id1]`; the bond style later forms `Omega_d = Omega - srot`.
  - **Y convention** (`subtract_groundstate=false`): `D = S^T R` (`Smat^T R`),
    `Phi_d = rotmat2euler(D)`, store in `fwd_euler[id1]`.
  - `fwd_Jinvtp[id1] = leftJacobianInverseTransposed(euler)`.
- Phase 4 forward-comms `fwd_euler` (3) + `fwd_Jinvtp` (9) to ghosts.

The **twist** is the third component:
`twist_wrapped_i = euler[2]` where `euler` is `Omega` (X) or `Phi_d` (Y).
Angle/dihedral styles read `fwd_euler[id1]`, `fwd_euler[id2]`, `fwd_euler[id3]`
(the left atom of each of their sub-steps) — so **everything is already keyed to the
left atom of the relevant step**, and a winding counter keyed the same way drops in
cleanly.

Define for step *i*:
- `phi_i`  = wrapped twist returned by the extractor, in `(-π, π]`.
- `n_i`    = integer winding counter (number of full ∓ turns accumulated).
- `Phi_i`  = **unwrapped twist** = `phi_i + 2π n_i`  (the quantity we want).

For the **Y convention** `phi_i` is the *deformation* twist (already small near the
duplex ground state because the intrinsic twist `S` is factored out). For the **X
convention** `phi_i = Omega_3` is the *absolute* step twist (≈ intrinsic ground-state
twist, e.g. ~0.6 rad/step at 1 bp, or ~6 rad ≈ near π folding at 10 bp). **The
melting model should use the Y convention**, because at 10 bp/bead the intrinsic
twist per bead is ~6 rad (>π) and would already sit *on the branch cut* of
`rotmat2euler` in the X convention; the Y convention removes that by factoring `S`
out at the group level (this is exactly why `subtract_groundstate=false` exists).

---

## A. Can the current representation trace twist past ±180°? Where exactly it breaks.

**No — not as stored.** `rotmat2euler` (`so3.h`, lines 98–155) returns
`theta = acos((tr R − 1)/2) ∈ [0, π]`, then `omega = (theta / (2 sinθ)) · vee(R−R^T)`.
The magnitude `|Omega| = theta` is **hard-capped at π**. A rotation matrix is an
element of SO(3) and **cannot store a winding number**: a 360° twist is literally the
identity matrix, indistinguishable from no twist. Concretely the failure modes are:

1. **Magnitude cap (the principal break).** Even decomposing into (tilt, roll, twist),
   the rotation *vector* magnitude is ≤ π. Once the true accumulated step twist passes
   π, `rotmat2euler` returns the supplementary/aliased value: as the real twist goes
   π → π+ε, the returned `Omega_3` jumps to ≈ −π+ε (a −2π discontinuity). A 10-bp
   melting event that untwists by ~2π therefore *wraps* and the stored `euler[2]` is
   useless as a continuous variable.

2. **θ ≈ π axis-extraction branch (lines 126–154).** When `sin θ → 0` at `θ ≈ π`,
   the general formula `coeff = θ/(2 sinθ)` blows up, so the code switches to a
   diagonal-based extraction `u_k = sqrt((R_kk+1)/2)` with the sign of the dominant
   component taken **positive** (`std::sqrt`), and the off-diagonal terms fixing the
   *relative* signs of the other two components. This branch returns an axis with an
   **undetermined overall sign** (`±u` give the same `R` at θ=π) — i.e. a genuine
   branch cut. In the X convention with 10-bp beads the *deformation-free* step can sit
   near θ=π, so the sign of the returned twist can flip frame-to-frame even without a
   physical event. In the Y convention `D = S^T R` is near identity (θ≈0, the
   *opposite* branch, lines 108–114, which is smooth), so the Y convention also *avoids
   the θ=π branch* during normal duplex fluctuations — another reason to use it.

3. **No half-turn distinction.** Because only `R` is stored, two states differing by an
   integer number of full turns are bit-identical. Continuity (winding) information
   exists only in the *trajectory*, never in a single frame. Therefore unwrapping
   **must** be done incrementally along the trajectory (and, where needed, along the
   chain), and the integer must be **persisted** between steps.

**Conclusion.** The rotation-vector/Y-convention representation is fine *for the
elastic force law* (which only ever needs the wrapped deformation; see §E) but is
fundamentally incapable of *storing* the winding number. We must add (i) a
continuity-safe twist *extractor* that does not rely on `rotmat2euler`'s capped
magnitude, and (ii) a persistent integer winding counter updated by trajectory
continuity.

---

## B. Robust continuity-based unwrapping scheme

### B.1 Continuity-safe wrapped-twist extractor

Replace naive reliance on `rotmat2euler[2]` with a dedicated *twist* extractor that
(a) is smooth through θ=π in the relevant matrix and (b) returns a value in `(−π, π]`
with a well-defined sign. Two equivalent realizations:

**(B.1a) Quaternion half-angle about the body twist axis (recommended, cheapest).**
The twist axis is the step's helical axis = the body z of the step frame. For a
rotation `D` (use `D = S^T R` in the Y convention, or `R` in X) with quaternion
`q = (w, x, y, z) = rotmat2quat(D)` (already in `so3.h`, lines 169–208), the
**swing–twist decomposition** about the z-axis gives the twist angle directly:

```
phi = 2 * atan2(z, w);            // twist about body z-axis, in (-2π, 2π]
phi = wrap_to_pi(phi);           // fold to (-π, π]
```

`atan2(z, w)` is smooth and sign-definite everywhere except the measure-zero set
`w = z = 0` (a pure 180° tilt/roll with zero twist — never reached in normal duplex
dynamics, and handled by the continuity check in B.2). This **does not use** the
capped `acos` magnitude and **does not** hit the θ=π diagonal branch of
`rotmat2euler`. Sign convention: positive `phi` = right-handed twist about +z of the
left frame (matches `Omega_3` for small angles).

Helper to add to `so3.h`:

```cpp
// Twist angle (about body z) from a rotation matrix, smooth & sign-definite,
// folded to (-pi, pi]. Uses the quaternion z,w components (swing-twist).
inline double twist_about_z(const double D[3][3]) noexcept {
    double q[4];
    rotmat2quat(D, q);            // q = (w, x, y, z), normalized
    double phi = 2.0 * std::atan2(q[3], q[0]);   // q[3]=z, q[0]=w
    // fold to (-pi, pi]
    const double TWO_PI = 2.0 * M_PI;
    while (phi >  M_PI) phi -= TWO_PI;
    while (phi <= -M_PI) phi += TWO_PI;
    return phi;
}
```

Note this twist agrees with `rotmat2euler`'s `omega[2]` to O(θ³) for small total
rotations (so the existing force law is unchanged when we keep using `euler[2]` for
forces — see §E), but it is robust at large tilt/roll where `rotmat2euler` mixes
components.

**(B.1b) Parallel-transport fallback (Bergou 2008 / Discrete Elastic Rods).** If the
swing–twist extractor proves fragile in tightly bent configurations, carry a Bishop
(zero-twist) reference frame transported bead-to-bead and define `phi_i` as the
signed angle between the transported material frame and the actual frame. This is the
mathematically cleanest unwrapped twist (additive by construction) but adds a framing
layer. **Recommendation:** ship B.1a as primary; keep B.1b documented as the fallback
toggled by a fix keyword.

### B.2 Per-junction winding update rule (trajectory continuity)

Store, per left atom *i*: `phi_prev_i` (wrapped twist last step) and `n_i` (integer
winding). Each `compute_lrf` call, after extracting the current wrapped `phi_i`:

```
dphi = phi_i - phi_prev_i;       // raw change since last frame
// detect a wrap: a physical step cannot change twist by ~2π in one MD step
if      (dphi >  PI) n_i -= 1;   // phi jumped UP by ~2π  => actually wound DOWN
else if (dphi < -PI) n_i += 1;   // phi jumped DOWN by ~2π => actually wound UP
phi_prev_i = phi_i;
Phi_i = phi_i + 2*PI*n_i;        // unwrapped twist (the CV / branch selector)
```

(Standard `np.unwrap` logic.) The sign mapping above follows from: if the *true*
unwrapped twist increased smoothly through +π, the extractor folds it to ≈ −π, i.e.
`dphi ≈ −2π`, so we **add** a turn — consistent with the comment. Verify the sign once
against a slow ramp test (Test T1).

**Sub-step safety (the one real correctness condition).** Trajectory unwrapping is
unambiguous iff `|true Δphi| < π` between consecutive *fix invocations*. The fix runs
every step in `pre_force`, so the relevant interval is **one MD step**, not the dump
interval — this is far safer than post-hoc trajectory unwrapping. At 10 bp/bead with
torsional relaxation times ≫ dt, a single step changes the deformation twist by
≪ 0.1 rad, so the condition holds with a huge margin. Two safeguards:

- **Guard band, not a hard threshold:** use `±π` exactly; the physics keeps |Δphi|
  tiny, so the choice of guard within (say) `[π/2, π]` is immaterial.
- **Optional jump alarm:** if `|dphi| > π_alarm` (e.g. 0.8π) *and* the configuration is
  not near `w=z=0`, emit a throttled warning (mirroring the FENE warning throttle in
  `bond_rbp_fene.cpp`, `last_fene_warn_step`), because that signals dt too large or an
  integrator blowup, both of which would also corrupt forces.

**What happens at the π branch cut.** With the swing–twist extractor (B.1a) there is
no θ=π diagonal branch involvement at all: `phi` from `atom2(z,w)` is continuous as the
*twist* passes π; the unwrap rule converts that continuity into an `n_i` increment.
The only genuine singularity is `w=z=0` (twist undefined because the step is a pure
π flip about an in-plane axis). At 10 bp this never occurs for the *deformation* `D`
in the Y convention (D ≈ identity in the duplex basin); in the melted basin the
deformation twist sits near `−2π·(per-bead turns)` but is *also tracked via n_i*, so
`D`'s wrapped part is again near identity. We therefore never sit on the singular set.

### B.3 What gets stored vs derived

- **Persisted (restart + exchange + ghost):** `phi_prev_i` (double), `n_i` (stored as
  double for buffer uniformity, holds an exact integer).
- **Derived each step (not persisted):** `Phi_i` (unwrapped twist), exposed read-only.

---

## C. Where to store it: extend `fix_rbp_lrf`

`fix_rbp_lrf` is the natural and only correct home: it already owns per-atom
persistent arrays with full `grow/copy/set/pack_exchange/forward_comm` plumbing, is
keyed to the left atom, and runs every step in `pre_force` (before forces). It
currently **lacks** `write_restart`/`read_restart` — we add them.

### C.1 New persistent per-atom arrays (`fix_rbp_lrf.h`)

```cpp
// existing:
double **triads;       // [nmax][9]
double **fwd_euler;    // [nmax][3]
double **fwd_Jinvtp;   // [nmax][9]

// NEW persistent winding state, keyed to the LEFT atom of each step:
double *tw_prev;       // [nmax]  wrapped twist last step, phi_prev_i
double *tw_wind;       // [nmax]  integer winding counter n_i (stored as double)
double *tw_unwrap;     // [nmax]  derived: Phi_i = phi_i + 2*PI*n_i (exposed)
int    tw_init_flag;   // 0 until the first compute_lrf seeds phi_prev/n
```

`tw_unwrap` is *derived* but kept in a per-atom array so a `compute` can read it and
it can be dumped; it need not be restarted (recomputed first step after restart from
`tw_prev`,`tw_wind`). `tw_prev` and `tw_wind` are the true persistent state.

### C.2 Lifecycle hooks (`fix_rbp_lrf.cpp`)

**Constructor:** allocate via `grow_arrays(atom->nmax)`; set `tw_init_flag = 0`. Note:
keep `comm->add_callback(Atom::GROW)` (already present). Add a **restart callback**:
`atom->add_callback(Atom::RESTART)` and set `restart_peratom = 1` so LAMMPS routes
per-atom restart through `pack_restart`/`unpack_restart` (see C.5).

**`grow_arrays(nmax)`** — add:
```cpp
memory->grow(tw_prev,   nmax, "fix_rbp_lrf:tw_prev");
memory->grow(tw_wind,   nmax, "fix_rbp_lrf:tw_wind");
memory->grow(tw_unwrap, nmax, "fix_rbp_lrf:tw_unwrap");
```

**`copy_arrays(i,j,delflag)`** — add: copy `tw_prev[j]=tw_prev[i]`, `tw_wind`, `tw_unwrap`.

**`set_arrays(i)`** (new atoms) — `tw_prev[i]=0; tw_wind[i]=0; tw_unwrap[i]=0;`
and mark "unseeded" so the first `compute_lrf` after creation seeds `phi_prev`
*without* incrementing `n` (see C.4). A clean way: store a sentinel `tw_prev[i]=NAN`
meaning "not yet seeded"; on the first pass detect NAN and seed instead of diff.

**`pack_exchange(i,buf)` / `unpack_exchange`** — append the three winding doubles.
Bump the count from 21 to **24**. (These migrate state when an atom changes proc.)

```cpp
// pack_exchange, after the existing 21:
buf[m++] = tw_prev[i];
buf[m++] = tw_wind[i];
buf[m++] = tw_unwrap[i];
return m;   // 24
```

### C.3 Forward comm to ghosts (`pack_forward_comm` / `unpack_forward_comm`)

Angle/dihedral styles read `fwd_euler[id1..id3]` on ghosts; a *compute* that sums
unwrapped twist also needs ghost values when the left atom is a ghost. **But the
winding update must happen exactly once per step on the owner of the left atom**, then
be *broadcast* to ghosts read-only. So:

- The **winding update (B.2) runs only for `id1 < nlocal`** (owned left atoms) inside
  Phase 2 of `compute_lrf`. Ghost left atoms are skipped for the update.
- Add `tw_unwrap` (and, if a downstream consumer needs them, `tw_prev`,`tw_wind`) to
  the **mode-0 forward comm**. Increase `comm_forward` from 12 to 13 (push
  `tw_unwrap` alongside `fwd_euler`+`fwd_Jinvtp`); ghosts then carry the correct,
  owner-computed unwrapped twist. `tw_prev`/`tw_wind` on ghosts are never used for the
  *update* (only owners update), so they need not be forward-communicated — but
  forward-comm of `tw_unwrap` is required so ghost steps report the right CV.

```cpp
// pack_forward_comm, comm_mode==0 branch, per atom j:
for (k<3)  buf[m++] = fwd_euler[j][k];
for (k<9)  buf[m++] = fwd_Jinvtp[j][k];
buf[m++] = tw_unwrap[j];     // NEW (comm_forward = 13)
```

**Ghost ordering subtlety.** `compute_lrf` Phase 2 iterates `neighbor->bondlist`,
which lists a bond only on the proc owning its **left** atom (`bondlist[bid][0]`).
Therefore every step's winding is updated on exactly one proc (the left-atom owner),
exactly once — no double counting, matching the existing keying contract documented in
`compute_lrf` ("each atom is the left endpoint of at most one RBP bond"). After the
update, Phase 4's forward comm propagates `tw_unwrap` to ghost copies. This is
identical to how `fwd_euler` is already handled, so the MPI consistency story is
unchanged.

### C.4 Initialization at t=0 / seeding (no spurious turn at startup)

On the **first** `compute_lrf` (or first time a given left atom is seen, detected by
the NAN sentinel / `tw_init_flag`):
```
phi_i = extract_twist(...);
if (tw_prev[id1] is unseeded) { tw_prev[id1] = phi_i; tw_wind[id1] = 0; }
tw_unwrap[id1] = phi_i + 2*PI*tw_wind[id1];
```
i.e. **seed `phi_prev = phi_now`, `n = 0`** so the first step produces `dphi=0` and no
winding increment. This defines the **zero of winding at simulation start** — the
unwrapped twist is measured *relative to the initial configuration*, which is the
physically meaningful reference for a melting (Δtwist) trajectory. If an absolute
linking number is required, add the known initial `Tw₀` offset in the *compute*
(§D), not here.

### C.5 Restart (`write_restart`/`read_restart` — currently absent)

Two distinct restart mechanisms exist in LAMMPS; the fix needs the **per-atom**
mechanism (the bond style's `write_restart` handles *coefficients*; the fix must
handle *per-atom winding state*).

Set in constructor: `restart_peratom = 1;` and `atom->add_callback(Atom::RESTART);`.
Implement:

```cpp
// number of per-atom doubles stored in restart (n_i, phi_prev_i) -> 2 (+1 length word)
int FixRBPLRF::pack_restart(int i, double *buf) {
  int m = 0;
  buf[m++] = 3;            // LAMMPS convention: first slot = # values incl. this
  buf[m++] = tw_prev[i];
  buf[m++] = tw_wind[i];
  return m;                // 3
}

void FixRBPLRF::unpack_restart(int nlocal, int nth) {
  double **extra = atom->extra;
  // skip to this atom's block
  int m = 0;
  for (int k = 0; k < nth; k++) m += static_cast<int>(extra[nlocal][m]);
  m++;                     // skip the count word
  tw_prev[nlocal]   = extra[nlocal][m++];
  tw_wind[nlocal]   = extra[nlocal][m++];
  tw_unwrap[nlocal] = tw_prev[nlocal] + 2.0*M_PI*tw_wind[nlocal];
}

int  FixRBPLRF::size_restart(int /*nlocal*/) { return 3; }
int  FixRBPLRF::maxsize_restart()            { return 3; }
```

(Exact signatures per the LAMMPS `Fix` per-atom restart API: `pack_restart`,
`unpack_restart`, `size_restart`, `maxsize_restart`, with `restart_peratom=1` and the
`Atom::RESTART` callback.) `tw_unwrap` is *derived* and need not be written — it is
reconstructed in `unpack_restart` and refreshed on the first post-restart
`compute_lrf`. We deliberately persist `tw_wind` (the integer) and `phi_prev` so that
a restart resumes the **same** unwrapped trajectory; without this, restarts would reset
winding to 0 and corrupt any running biased simulation or hysteresis measurement.

The bond/dihedral/angle styles already restart their *coefficients*
(`BondRBP::write_restart`, `BondRBPFene::write_restart`); we do not touch those.

### C.6 `memory_usage`

Add `3 * nmax * sizeof(double)` for the three new arrays.

---

## D. Exposing unwrapped twist as a collective variable / compute

Add a thin **`compute rbp/twist`** that reads the fix's `tw_unwrap` array. Two outputs:

1. **Per-junction unwrapped twist** (a per-atom vector, keyed to the left atom):
   `compute ID group rbp/twist` → `c_ID[i]` = `Phi_i` for the step whose left atom is i
   (0 for atoms that are not a left endpoint of an RBP bond).
2. **Global accumulated twist / linking proxy** (a scalar):
   `Phi_total = Σ_i Phi_i` over all owned left atoms, MPI-`Allreduce`d. This is the
   total relative twist `Tw − Tw₀` of the chain (a writhe-free linking-number proxy in
   the *straight/clamped* reference; for the full Lk add a writhe term, below).

Implementation sketch (`compute_rbp_twist.cpp`):

```cpp
void ComputeRBPTwist::compute_vector() {       // global: [0]=Tw_total, [1]=<Phi>, ...
  FixRBPLRF *fix = /* modify->get_fix_by_style("^rbp/lrf")[0] */;
  double local = 0.0; int **bondlist = neighbor->bondlist;
  for (bid) { int id1 = bondlist[bid][0];
              if (id1 < atom->nlocal) local += fix->tw_unwrap[id1]; }
  MPI_Allreduce(&local, &vector[0], 1, MPI_DOUBLE, MPI_SUM, world);
}
void ComputeRBPTwist::compute_peratom() {       // per-atom Phi_i
  for (i<nlocal) vector_atom[i] = is_left_atom(i) ? fix->tw_unwrap[i] : 0.0;
}
```

(Iterate the **bondlist** rather than all atoms so we count each step once and respect
the left-atom keying; the fix has already filled `tw_unwrap[id1]` for owned and ghost
left atoms.)

**Use for biasing.** The global scalar `c_ID[0]` (`Tw_total`, unwrapped) is the natural
1-D reaction coordinate for torque-melting. Feed it to:
- `fix addtorque` / a custom torque term as the conjugate variable, or
- **PLUMED** via the LAMMPS-exposed compute, or a PLUMED `CUSTOM`/`MATHEVAL` CV — note
  PLUMED's built-in `TORSION` is wrapped to (−π,π] and **must not** be used for
  multi-turn twist; bias the unwrapped compute instead.
- Alternatively bias a **periodicity-free discrete CV** (count of "melted" steps,
  defined by `n_i ≤ −1` or `Phi_i < threshold`), which sidesteps angle periodicity
  entirely (recommended for FFS / committor analysis, per the review's
  enhanced-sampling section).

**Linking-number cross-check (optional, analysis only).** Provide a separate
`compute rbp/writhe` (Klenin–Langowski exact polygonal writhe over bead centerlines,
O(N²); or Fuller single-integral vs a straight reference for speed, falling back to
the double sum near antipodal points per Neukirch–Starostin 2008). Then
`Lk = Tw_total/2π + Wr` should be conserved at fixed linking number — a strong
per-frame validation of the whole twist-tracking machinery (White 1969).

---

## E. Force / torque consistency — adding +2πn shifts only the branch, not the Jacobian

**Claim.** Tracking `n_i` and forming `Phi_i = phi_i + 2π n_i` has **zero effect** on
the existing elastic forces/torques, and only changes *which branch of a non-quadratic
twist potential applies*. Justification, pinned to the code:

1. **The harmonic elastic law never sees `n`.** In `bond_rbp.cpp` / `bond_rbp_fene.cpp`
   the force/torque are built from `Omd` (= `fwd_euler[id1]`, the *wrapped* deformation)
   and the Jacobian `fwd_Jinvtp[id1] = leftJacobianInverseTransposed(euler)`. Adding
   `2πn` to a *reported* twist does not alter `R`, `D`, `euler`, or `Jinvtp` — those are
   functions of the rotation matrices only. The harmonic potential is, and stays, a
   function of the wrapped deformation. **We do not feed `Phi_i` back into the harmonic
   term.** `n_i` is bookkeeping layered *on top of* the unchanged `R`. So the existing
   2-pt/3-pt/4-pt couplings (off-diagonal `M` blocks in bond/angle/dihedral) are
   untouched.

2. **The Jacobian is a property of the manifold, not the chart winding.** The map from
   se(3) gradient to body torque, `torque = T1 · Jinvtp · (∂E/∂Omega)` (and the Y-conv
   `Smat`-sandwiched variant), depends only on the *local* rotation vector `Omega`
   (magnitude ≤ π by construction of the chart at this frame). `leftJacobianInverseTransposed`
   is evaluated at the wrapped `euler`, which is exactly the correct tangent-space
   representative. Shifting the *reported* twist by `2πn` does not change the tangent
   space or `Jinvtp`. (At θ→π the chart is singular, but as argued in §A/B the Y
   convention keeps `D` near identity, far from that singularity.)

3. **A double-well / scalar twist potential DOES use `Phi_i` — only to pick the branch.**
   When we later add the bimodal melted potential (a one-sided scalar twist potential,
   architecturally the `bond_rbp_fene` template: `E(twist)` on top of the wrench, with
   analytic `dE/dtwist` mapped to torque via the *same* `Jinvtp`), the *unwrapped*
   `Phi_i` selects whether we are in the duplex basin (`n_i = 0`) or the melted basin
   (`n_i = −1`, shifted by −2π). The **gradient** `dE/dΦ = dE/dφ` (since `Φ = φ + 2πn`
   with `n` piecewise-constant) — i.e. **the torque is computed from the same local
   `∂E/∂twist` and the same `Jinvtp`**; the winding only changes *which well's
   center/stiffness* enters. Thus `+2πn` shifts the energy branch but the
   force-mapping algebra is byte-identical to the FENE add-on. This is the precise sense
   in which winding is "free" for forces.

**Sign conventions to pin down (do once, lock with a unit test):**
- `phi_i > 0` ⇒ right-handed (positive) twist about +z of the **left** step frame,
  matching `Omega_3` for small angles and the existing `srot[2]`/`Ystatic[2]` twist
  component sign. Underwinding (melting) ⇒ `phi` decreasing ⇒ eventually `n_i = −1`.
- `dphi < −π ⇒ n += 1`, `dphi > +π ⇒ n -= 1` (B.2). Confirm against the duplex
  intrinsic right-handedness so that a +turn applied torque raises `Tw_total`.
- Applied torque term `−τ·Tw_total` with `τ < 0` (underwinding, ≈ −10 pN·nm melting
  torque) lowers the melted (`n=−1`) branch — verify the global scalar increases under
  positive applied torque.

---

## F. Edge cases

**F1. Closed (circular) topology** (`parse_rbp.h` metadata `closed`). The chain has a
"wrap-around" step closing the ring; its left atom is the last bead, bonded to the
first. Nothing special is needed for *trajectory* unwrapping (each step is still keyed
to its left atom and updated independently). But for the **global CV**, on a closed
loop the *total* `Σ phi_i` (wrapped) is the integer-multiple-of-2π topological twist,
and `Σ Phi_i` (unwrapped, relative to t=0) is the change in Tw. The linking number
`Lk = Tw + Wr` is then a true topological invariant — use the `rbp/writhe` compute to
verify `Lk` is conserved at fixed-Lk on the circle (the cleanest validation; White
1969, Klenin–Langowski 2000). Ensure the closing step is included exactly once in the
bondlist (it is, as an ordinary bond) and that the minimum-image handling in the bond
styles (`domain->minimum_image`) is already applied (it is, see `bond_rbp.cpp`
lines 144–146).

**F2. Initialization at t=0** — see §C.4. Seed `phi_prev=phi_now, n=0`; winding is
*relative* to the start. For an absolute Lk, add a constant `Tw₀` (known from the build)
in the compute, never in the per-atom state.

**F3. Parallel domain decomposition.** Covered in §C.3: the update runs only on the
owner of the left atom (`id1 < nlocal`), exactly once per step (bondlist lists each
bond once, on the left-atom owner), then `tw_unwrap` is forward-comm'd to ghosts.
`pack_exchange`/`unpack_exchange` carry `phi_prev`,`n`,`unwrap` when an atom migrates,
so winding survives load balancing. **Critical ordering guarantee:** between two
`compute_lrf` calls an atom may migrate procs, but its `phi_prev` migrates with it
(pack_exchange), so the next-step `dphi = phi_now − phi_prev` is still computed against
the correct previous value regardless of which proc now owns it. This is why
`phi_prev` *must* be in `pack_exchange` (not just forward-comm). Reneighboring/exchange
happens between steps, before `pre_force`, so the migrated `phi_prev` is in place.

**F4. Multiple bonds per atom.** The existing keying contract assumes **each atom is the
left endpoint of at most one RBP bond** (documented in `compute_lrf`, lines 498–500).
The winding arrays inherit this assumption — `tw_prev[id1]`/`tw_wind[id1]` belong to
the unique step whose left atom is `id1`. For a linear or simple circular backbone this
holds. If a future topology violates it (branch points), the design must switch from
"keyed to left atom" to "keyed to a per-bond index" (a parallel array indexed by
`bid` with its own exchange/restart) — flag this as out-of-scope for now and add a
hard check in `validate_junctions` (count RBP bonds per left atom; error if >1).

**F5. Interaction with 3-pt / 4-pt couplings.** The angle/dihedral styles read the
*wrapped* `fwd_euler[id1..id3]` for forces (unchanged), and only the *new compute* and
the *future double-well term* read `tw_unwrap`. Because §E guarantees winding does not
touch the harmonic/coupled force law, the banded 2/3/4-pt couplings are unaffected.
The double-well term (when added) is a per-step scalar twist potential (FENE-template),
local to one step; it does not enter the off-diagonal couplings, so no
cross-term bookkeeping is needed beyond what the existing styles already do.

**F6. Energy minimization (`min_pre_force`).** `compute_lrf` also runs during
minimization. During a minimize the "trajectory" is a sequence of trial geometries
that can jump by more than π (line searches), which would *spuriously* increment `n`.
**Guard:** disable winding updates when `update->whichflag == 2` (minimization) — i.e.
in `min_pre_force` call a variant that refreshes `phi_prev` and `tw_unwrap` but does
**not** diff/increment `n`. Winding tracking is only meaningful along an MD trajectory.

---

## G. Implementation plan (step by step) and validation

### G.1 Code changes (ordered, minimal-risk)

1. **`so3.h`:** add `twist_about_z(const double D[3][3])` (B.1a). Pure addition, no
   change to existing functions. (Optionally add `wrap_to_pi`.)
2. **`fix_rbp_lrf.h`:** declare `tw_prev, tw_wind, tw_unwrap` (public, so the compute
   and future double-well term can read them, mirroring `fwd_euler`); declare
   `pack_restart/unpack_restart/size_restart/maxsize_restart`; bump `comm_forward` note.
3. **`fix_rbp_lrf.cpp`:**
   - Constructor: `comm_forward = 13;` `restart_peratom = 1;`
     `atom->add_callback(Atom::RESTART);` allocate new arrays; init to NAN-sentinel.
   - `grow_arrays/copy_arrays/set_arrays/memory_usage`: extend for 3 arrays.
   - `pack_exchange/unpack_exchange`: 21 → 24.
   - `pack_forward_comm/unpack_forward_comm` (mode 0): append `tw_unwrap`.
   - `compute_lrf` Phase 2: after computing `euler`, compute
     `phi = twist_about_z(D_or_R)`; if `id1 < nlocal` run the seed-or-update rule
     (B.2/C.4) and set `tw_unwrap[id1]`. Skip update for ghost left atoms.
   - Add `pack_restart/unpack_restart/size_restart/maxsize_restart` (C.5).
   - Add a `min_pre_force` path that seeds-without-incrementing (F6).
   - Destructor: `atom->delete_callback(id, Atom::RESTART);` `memory->destroy` new arrays.
4. **`compute_rbp_twist.{h,cpp}`** (new): per-atom `Phi_i` + global `Tw_total` (D).
5. *(Later, separate task)* **`compute_rbp_writhe`** and the **double-well scalar twist
   potential** (FENE-template) consuming `tw_unwrap` as the branch selector.
6. **`install.sh`:** add the new compute source(s).

No changes to `bond_rbp.cpp`, `dihedral_rbp.cpp`, `angle_rbp.cpp` forces — winding is
purely additive (§E).

### G.2 Validation tests

- **T0 (extractor unit test).** For random `D ∈ SO(3)`, check `twist_about_z(D)` equals
  `2·atan2(qz,qw)` folded, and agrees with `rotmat2euler(D)[2]` to O(θ³) for small θ.
  Check smoothness sweeping twist from 0→2π in a *fixed* tilt/roll: `phi` is continuous
  and the unwrap rule yields a clean ramp `Phi: 0→2π` with one `n` increment.
- **T1 (single-step slow ramp, sign lock).** One step, apply a constant torque to one
  end (`fix addtorque`) and clamp the other; ramp `Tw_total` past +1 and −1 turn at
  small dt. Assert `Tw_total` is monotone and continuous through ±π (no sawtooth), and
  fix the sign convention (positive applied torque ⇒ `Tw_total` increases ⇒ `n` rises).
- **T2 (Lk conservation, closed circle).** Build a relaxed minicircle (`closed`),
  `fix rigid`/clamped, fixed Lk. Each frame compute `Tw_total/2π` (this design) and
  `Wr` (Klenin–Langowski). Assert `Tw + Wr = Lk = const` to numerical tolerance across
  a buckling event — the decisive end-to-end test (White 1969).
- **T3 (MPI invariance).** Same trajectory on 1 vs N procs with aggressive
  decomposition and load balancing; assert per-junction `Phi_i` and global `Tw_total`
  agree bit-for-bit-modulo-roundoff across decompositions and across a migration event
  (atom crossing a proc boundary mid-trajectory).
- **T4 (restart round-trip).** Run M steps, write restart, continue; vs run 2M steps
  uninterrupted. Assert `Tw_total(t)` and all `n_i`, `phi_prev_i` match after restart
  (proves C.5). Specifically check a restart taken *mid-winding* (just after an `n`
  increment) resumes without resetting winding.
- **T5 (force non-perturbation).** Confirm total force, torque, and potential energy are
  **bit-identical** with vs without the winding code compiled in (a no-op for the
  harmonic law; §E). This guards against accidental coupling.
- **T6 (melting event).** Under underwinding torque past the melting plateau, verify a
  ~1-turn untwist registers as `n_i: 0 → −1` on the melted step(s) and `Tw_total`
  drops by ~2π, while `Lk` (T2 check) stays fixed because writhe compensates — the
  physical signature the whole machinery exists to capture (Sheinin–Wang 2011;
  Matek et al. 2015).
- **T7 (minimize guard).** Run `minimize`; assert no `n` increments occur during line
  searches (F6), then MD resumes tracking correctly.

### G.3 Risks / mitigations

- *Spurious winding from large per-step jumps:* mitigated by per-MD-step (not
  per-dump) unwrapping + jump alarm (B.2). If ever triggered legitimately, reduce dt.
- *θ=π singularity:* avoided by Y convention (`D≈I`) + swing–twist extractor (§A/B).
- *Multiple-bond-per-atom topologies:* hard-error guard in `validate_junctions` (F4).
- *Restart API mismatch:* the per-atom restart path (`restart_peratom`+`Atom::RESTART`)
  is the only one that migrates per-atom winding; verify signatures against the LAMMPS
  version in tree before merge.

---

## Summary

The rotation-matrix / Y-convention representation **cannot store winding** (|Omega|≤π
cap; θ=π branch), but it is fine for forces, which only need the wrapped deformation.
Add a continuity-safe **swing–twist extractor** (`twist_about_z`, quaternion
half-angle — no `acos` cap, no θ=π branch), a per-left-atom **winding counter `n_i`**
with `phi_prev_i`, updated **once per MD step on the left-atom owner** by the standard
±2π-jump rule, seeded at t=0 to zero winding, and propagated to ghosts via the existing
forward comm. Store the state in **`fix_rbp_lrf`** (extend
grow/copy/set/exchange/forward arrays; **add the missing per-atom
restart**), expose it via a new **`compute rbp/twist`** (per-junction `Phi_i` + global
`Tw_total`) for biasing/analysis, and cross-check topology with `Lk = Tw + Wr`. Adding
`+2πn` shifts only **which branch of a (future) double-well twist potential applies**;
the Jacobian/torque mapping and the existing 2/3/4-pt harmonic couplings are provably
unchanged. Validation centers on Lk conservation on a closed minicircle, MPI/restart
invariance, and bit-identical forces.
