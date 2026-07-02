# Phenomenological Bimodal (Duplex / Melted) Free-Energy Functional for the SE(3) CG-RBP Model

**Design document — torque-induced melting as barrier crossing at ~10 bp / bead.**

This document specifies a concrete, implementation-ready bimodal free-energy term to bolt onto the
existing CG-RBP LAMMPS force field. It is grounded in the actual source
(`bond_rbp_fene.cpp`, `so3.h`, `fix_rbp_lrf.cpp/.h`, `parse_rbp.h`) and in the literature synthesis
(Sicard–Manghi 2015; Marko–Neukirch 2013; Marko 2007; Sheinin–Wang 2011; Cocco–Monasson–Marko 1999;
SantaLucia 1998; Wiggins–Phillips–Nelson 2005; Argudo–Purohit 2014; Dans–Orozco 2020/2021;
ICNN log-sum-exp arXiv:2506.17242).

---

## 0. Notation and the existing harmonic baseline

Per RBP step *i* (a bond between left bead `id1` and right bead `id2`), the existing model builds a
6-vector deformation `Yd` from the relative junction `g_i = T_i^T T_{i+1}`:

```
Yd = (Phi_d, w_d)            Phi_d in R^3 (rotational), w_d in R^3 (translational)
```

- **X convention** (`subtract_groundstate = true`): `Phi_d = Om - srot`, where `Om = log(T1^T T2)`
  (`so3::rotmat2euler`), and `w_d = w - svec`, `w = T1^T (r2 - r1)`.
- **Y convention** (`subtract_groundstate = false`): `Phi_d = log(S^T R)` with `S = euler2rotmat(srot)`
  precomputed (`Smat`), and `w_d = S^T (w - svec)`.

In both conventions `Yd` is a *deformation away from the ground state* and is *small* in the duplex
basin. The energy is the single harmonic well

```
beta*E_dup = (1/2) Yd^T M Yd                                              (Eq. 1)
```

with `M` the 6x6 stiffness (banded across steps via the 2-pt/3-pt/4-pt couplings in the bond/angle/
dihedral styles). The third rotational component `Phi_d[2] = Yd[2]` is the **twist deformation**;
`Omega_3 = log(R)_3` is the bare twist, capped at pi by `rotmat2euler`. Forces/torques are obtained
by mapping the se(3) gradient `dE/dYd` to body torques with `so3::leftJacobianInverseTransposed`
(see `bond_rbp_fene.cpp` lines 171-253). This is the architecture the bimodal term must mirror.

**Key precomputed quantities available from `fix_rbp_lrf`** (keyed to the left atom `id1`):
`triads[i]`, `fwd_euler[i]` (= `Om` for X, `Phi_d` for Y), `fwd_Jinvtp[i]`
(= `leftJacobianInverseTransposed` of that euler vector). The bimodal style reads these via zero-copy
pointer cast exactly as `BondRBPFene::compute` already does.

**Challenge-1 dependency (the unwrapped twist).** The melted basin is shifted by ~ -1 turn in twist;
a basin centered at `-2pi` cannot be represented by a rotation-vector component that lives in
`[-pi, pi]`. We therefore require a **continuous, unwrapped per-step twist** `phi_unw(i)` (winding
number tracked across the trajectory). This is to be produced by the Challenge-1 extension of
`fix_rbp_lrf` (per-step winding counter + restart). The bimodal term *consumes* `phi_unw(i)`; its
production is specified in §7.2 but is out of scope here.

Throughout, energies are in units of `kBT` (the model already works in reduced/beta units, cf. the
`beta*E` energy tallies). `beta = 1`.

---

## 1. Candidate functional forms

We compare four families. Let

```
E_dup(Yd)  = (1/2) (Yd - Yd0_dup )^T M_dup  (Yd - Yd0_dup )          (stiff duplex well)
E_melt(Yd) = (1/2) (Yd - Yd0_melt)^T M_melt (Yd - Yd0_melt) + dG     (soft melted well + offset)
```

In the **Y convention the duplex ground state is the origin**, `Yd0_dup = 0`, because the static `S`
is already factored out. The melted ground state differs from duplex *only* in the components we
choose to shift; the dominant shift is in the unwrapped twist:

```
Yd0_melt = Yd0_dup + Delta,   Delta = (0, 0, dphi_melt, 0, 0, drise_melt)               (Eq. 2)
```

with `dphi_melt ≈ -2pi` per 10-bp bead (one left-handed turn; Sheinin–Wang 2011 L-DNA ≈ -13 bp/turn)
and optional rise extension `drise_melt` (L-DNA rise 0.46-0.48 nm/bp vs 0.34). `M_melt << M_dup`
(soft bend/twist, §5). **`dphi_melt` lives in the unwrapped-twist coordinate, not the wrapped
rotation vector** — this is the crux of why Challenge 1 is a prerequisite.

### (i) Log-sum-exp / soft-min of two wells (RECOMMENDED)

```
beta*E(Yd) = -(1/rho) * ln[ exp(-rho * E_dup) + exp(-rho * E_melt) ]                     (Eq. 3)
```

- Smooth, differentiable everywhere; one extra scalar `rho` controls how sharply the two wells are
  glued. `rho -> inf` => hard `min(E_dup, E_melt)` (cusp + cusped barrier ridge);
  `rho` small => merged single well.
- The special case `rho = beta = 1` is exactly the **two-state Boltzmann sum** of Wiggins–Phillips–
  Nelson (KWLC): `E = -kT ln[exp(-beta E_dup) + exp(-beta E_melt)]`. That is the physically motivated
  default (it is a genuine marginal free energy over a hidden 2-state spin), and it is also the k=2
  Gaussian-mixture of Dans–Orozco 2020 with weight ratio `exp(-beta*dG)`.
- **Barrier height is NOT independent of the wells** at fixed `rho`: the saddle sits at the
  intersection `E_dup = E_melt` and is rounded down by `~ (1/rho) ln 2`. To make the kinetic barrier
  an **independent knob**, add an explicit Gaussian bump on the unwrapped twist at the crossover
  (§4).
- Gradient is a softmax-weighted average of the two well gradients — a few extra lines on top of the
  existing harmonic gradient (§3). **This is the recommended primary form.**

### (ii) Discrete per-step Ising spin `s_i in {0,1}` (Poland–Scheraga / Benham / Storm–Nelson)

```
H = sum_i [ (1 - s_i) E_dup,i + s_i (E_melt,i + dG_i(seq)) ]
        + J * sum_i (s_i != s_{i+1})            (domain-wall / nucleation cost)
        - tau * Tw                              (applied torque tilt)                    (Eq. 4)
```

- Cleanest encoding of cooperativity: the junction penalty `J` is literally the nucleation barrier;
  equilibrium occupancy and the constant-torque plateau follow from a 2x2 transfer matrix
  (Storm–Nelson 2003; Benham 1992; Fye–Benham 1999 for the global Lk constraint).
- Forces in MD require either (a) a hybrid **MC move on `s_i`** interleaved with MD (Metropolis on the
  per-step energy difference, cheap because `s_i` is local), or (b) promotion of `s_i` to a continuous
  auxiliary `lambda_i in [0,1]` with its own (over-damped) dynamics — *lambda-dynamics*:
  `E(Yd, lambda) = (1-h(lambda)) E_dup + h(lambda)(E_melt + dG) + V_bias(lambda)`, switching function
  `h(lambda) = lambda^2(3-2lambda)` (smooth, `h'=0` at endpoints), and a double-well biasing
  potential `V_bias` on `lambda` to keep it bimodal. This is the Gaussian-mixture (iv)/(i) in
  disguise once `lambda` is integrated out.
- **Verdict:** keep as the *cooperativity layer* (the `J` term, mapped onto the existing 3-pt/4-pt
  couplings, §6), but not as the per-step MD force law, because a discrete spin needs an extra
  integrator/MC machinery the current pure-MD pipeline lacks.

### (iii) Landau quartic in the unwrapped twist

```
beta*V(phi) = A4 (phi - phi0)^4 + A2 (phi - phi0)^2 + A1 (phi - phi0),   A2 < 0,  phi = Tw_unw   (Eq. 5)
```

- Minimal analytic double well: minima at `phi0 ± sqrt(-A2/2A4)`, barrier `A2^2/(4 A4)`, linear tilt
  `A1 = -tau * (dTw/dphi)` from applied torque (Argudo–Purohit 2014; Landau). Ideal for the reduced
  1-D Kramers/PMF analysis.
- **Limitation:** the two well curvatures are locked together by `(A2, A4)`, so it *cannot*
  independently set a stiff duplex and a soft melted curvature, and it acts on twist alone (no
  bend-softening). Use it for the reduced-coordinate validation (§8), not as the full 6-D step energy.

### (iv) Other literature forms considered

- **Per-phase parabolas + Maxwell construction** (Marko 2007; Marko–Neukirch 2013): not an MD force
  law but the **calibration constraint** (§5.4): the coexistence torque
  `tau_c = (dG)/(phi_dup - phi_melt)` must equal -10..-11 pN·nm.
- **Morse + anharmonic stacking (PBD)**: introduces an extra opening DOF per bead; unnatural at 10 bp.
  We borrow only its central lesson — *the melted state is a softened-stiffness state, not merely an
  energy offset* (Dauxois–Peyrard–Bishop 1993).
- **Sicard–Manghi 2015 distance-dependent torsional modulus** `kappa_phi(rho): 450 kBT -> 0`: this is
  exactly the `M_dup -> M_melt` switch, and their result that the **barrier scales affinely with
  beta·kappa_phi** is our justification for treating `M_melt` (and the bump `B`) as the tunable
  barrier knobs.

**Recommendation:** **(i) log-sum-exp of two shifted Gaussian wells** as the per-step MD energy
(with `rho = 1` default), **plus an explicit twist barrier bump** (§4) for independent kinetics,
**plus the Ising junction penalty (ii) realized through the existing 3-pt/4-pt M couplings** (§6) for
cooperativity, calibrated by the Maxwell construction (iv) and reduced to the quartic (iii) for
Kramers analysis.

---

## 2. Recommended energy — full expression

Per step *i*, with unwrapped twist `phi = phi_unw(i)` substituted into the rotational-twist slot of
`Yd` (see §2.1):

```
beta*E_i(Yd)  =  E_LSE  +  E_bump                                                        (Eq. 6)

E_LSE  = -(1/rho) ln[ exp(-rho * E_dup) + exp(-rho * E_melt) ]

E_dup  = (1/2) Yd^T M_dup Yd                                  (Yd0_dup = 0, Y convention)
E_melt = (1/2) (Yd - Delta)^T M_melt (Yd - Delta) + dG(seq)

E_bump = B * exp( -(phi - phi_b)^2 / (2 w_b^2) )             (kinetic barrier, twist-only)
```

### 2.1 Mixing the unwrapped twist into `Yd`

`Yd` as computed by the existing code has its twist component `Yd[2]` wrapped to `(-pi, pi]`. Define
the **augmented deformation** `Yd*` identical to `Yd` except in the twist slot:

```
Yd*[k] = Yd[k]            for k != 2
Yd*[2] = phi_unw(i) - srot_twist_dup        (X conv: subtract duplex intrinsic twist)
Yd*[2] = phi_unw(i)                         (Y conv: S already removed, so phi_unw is the deformation)
```

`phi_unw(i)` is the winding-aware twist from `fix_rbp_lrf` (Challenge 1). In the duplex basin
`phi_unw` differs from the wrapped `Yd[2]` by a multiple of `2pi` that is **zero** (no winding), so
`Yd* == Yd` exactly there — the harmonic baseline is recovered (§2.3). Only once a step has wound by
~ -1 turn does `Yd*[2] ≈ -2pi` while the wrapped `Yd[2]` would have jumped — `Yd*` is the physically
correct, continuous coordinate. All `E_dup`, `E_melt`, `E_bump` below use `Yd*`; we drop the star.

### 2.2 Parameters and physical meaning

| symbol | meaning | typical value / source |
|---|---|---|
| `M_dup` (6x6) | duplex stiffness | existing `Mmat` (cgNA+/PolyCG at 10 bp; A≈45-50 nm, C≈95-110 nm) |
| `Yd0_dup` | duplex ground state | `0` in Y convention (S factored out) |
| `M_melt` (6x6) | melted (L-DNA / ssDNA) stiffness | soft: bend persistence ≈3-7 nm, twist modulus ≈19-20 nm (Sheinin–Wang 2011; Marko–Neukirch 2013), or ssDNA-limit ≈1-3 nm / ≈0 |
| `Delta` (6-vec) | duplex->melted ground-state shift | `Delta[2] = dphi_melt ≈ -2pi` per bead; `Delta[5] = drise_melt ≈ +0.12-0.14 nm/bp * 10` optional |
| `dG(seq)` | sequence-dependent duplex stabilization (offset between basin minima) | `+7..+24 kBT` per 10-bp bead from SantaLucia NN (§5.3) |
| `rho` | soft-min sharpness | default `1` (= KWLC Boltzmann sum); raise to sharpen crossover |
| `B` | barrier bump height | tune to ~22 kBT nucleation / ~13 kBT closure (Sicard–Manghi 2015); see §4 |
| `phi_b` | bump center (twist) | midway between basins, `≈ -pi` (i.e. `dphi_melt/2`) |
| `w_b` | bump width (twist) | `≈ 0.5-1.0 rad`; sets curvature of the transition state |
| `tau` | applied torque (optional Legendre tilt) | added as `-tau * phi`, §7.4 |

### 2.3 Exact reduction to the existing harmonic well

In the duplex basin `Yd` is small and `phi_unw = 0` (no winding), so `Yd* = Yd`:

```
E_melt - E_dup = (1/2)(Yd-Delta)^T M_melt (Yd-Delta) - (1/2) Yd^T M_dup Yd + dG
```

With `|Delta[2]| = 2pi` and `dG > 0`, evaluated at `Yd ≈ 0` this is `≈ (1/2) Delta^T M_melt Delta + dG`,
a large positive number (tens of `kBT`). Hence `exp(-rho E_melt) << exp(-rho E_dup)` and

```
E_LSE = E_dup - (1/rho) ln[ 1 + exp(-rho (E_melt - E_dup)) ]
      = E_dup - (1/rho) * exp(-rho (E_melt - E_dup)) + O(exp(-2 rho dE))
      -> E_dup = (1/2) Yd^T M_dup Yd                                                     (Eq. 7)
```

with exponentially small correction. Setting `B = 0` (no bump, or `phi_b` far from the duplex basin so
`E_bump ≈ 0` there) gives **bit-for-bit the existing energy in the duplex basin**, and the gradient
likewise reduces to `M_dup Yd` (§3). This guarantees the new style is a strict superset: with
`dG -> +inf` (or melted basin disabled) it *is* `bond_rbp` / `bond_rbpfene`. This is the backward-
compatibility contract for validation.

---

## 3. Analytic gradient and force/torque mapping

Let `g_dup = M_dup * Yd` and `g_melt = M_melt * (Yd - Delta)` be the two harmonic gradients
(`dE/dYd`, each a 6-vector). Define softmax weights

```
p_dup  = exp(-rho E_dup ) / Z,   p_melt = exp(-rho E_melt) / Z,   Z = exp(-rho E_dup)+exp(-rho E_melt)
```

Then the gradient of the log-sum-exp is the **convex combination**

```
dE_LSE/dYd = p_dup * g_dup + p_melt * g_melt                                             (Eq. 8)
```

(The `rho` factors cancel exactly — this is the standard softmax-of-quadratics gradient. Note
`p_dup + p_melt = 1`, so in the duplex basin `dE_LSE/dYd -> g_dup = M_dup Yd`, matching Eq. 1.)

The bump contributes only to the twist component:

```
dE_bump/dYd[2] = -B (phi - phi_b)/w_b^2 * exp(-(phi-phi_b)^2/(2 w_b^2))                  (Eq. 9)
dE_bump/dYd[k] = 0  (k != 2)
```

So the **total se(3) gradient** is

```
A_total = (dE/dPhi_d) = rotational part (indices 0,1,2) of (p_dup g_dup + p_melt g_melt) + [0,0, dE_bump]
B_total = (dE/dw_d)   = translational part (indices 3,4,5) of (p_dup g_dup + p_melt g_melt)
```

**Force/torque mapping is identical to `bond_rbp_fene.cpp`** — only the (A, B) vectors change. Reusing
the Y-convention block (lines 223-253), with `Smat` and `Jinvtp` read from `fix_rbp_lrf->fwd_Jinvtp`:

```
torque_2 = T1 * Smat * (Jinvtp * A_total)
force_1  = T1 * Smat * B_total
torque_1 = T1 * (w x (Smat * B_total)) + torque_2
```

and the equal-and-opposite application to `id1`/`id2` (lines 370-391). **Crucial correctness point:**
`Jinvtp` in the current code is `leftJacobianInverseTransposed(Phi_d)` evaluated at the *wrapped*
euler vector. The left Jacobian is `2pi`-periodic in the rotation angle in the sense that the
*physical* body torque depends only on the actual orientation, which `fwd_euler` (wrapped) already
encodes correctly. Because `A_total`'s twist component is the gradient *value* (a scalar magnitude,
finite even when `phi` is large), and the orientation enters only through the physical triads, the
mapping is valid: **the unwrapped twist enters the energy/gradient magnitude, but the geometric
torque-axis mapping uses the true (wrapped) orientation.** The bump derivative `dE_bump/dphi` is a
torque about the local twist axis and is mapped through the same `Jinvtp` row.

Pseudo-code for the compute kernel (drop-in replacement for the wrench block):

```cpp
// --- read precomputed step data (zero-copy), Y convention ---
double (*T1)[3]    = (double(*)[3]) fix_lrf->triads[id1];
double Phi_d[3]    = { fwd_euler[id1][0], fwd_euler[id1][1], fwd_euler[id1][2] }; // wrapped
double (*Jinvtp)[3]= (double(*)[3]) fix_lrf->fwd_Jinvtp[id1];
double (*Smat)[3]  = params[bt].Smat;

// --- translational deformation (unchanged) ---
double w[3]; lamath::mul_Atx(T1, dr, w);
double tmp[3]; lamath::subtract(w, params[bt].svec, tmp);
double wd[3];  lamath::mul_Atx(Smat, tmp, wd);

// --- augmented Yd with UNWRAPPED twist (Challenge 1) ---
double phi_unw = fix_lrf->twist_unwrapped[id1];   // continuous, winding-aware
double Yd[6] = { Phi_d[0], Phi_d[1], phi_unw, wd[0], wd[1], wd[2] };
//                              ^^^^^^^ replaces wrapped Phi_d[2]

// --- two harmonic wells ---
double Ydm[6];  for (k=0;k<6;k++) Ydm[k] = Yd[k] - params[bt].Delta[k];
double Edup  = 0.5 * quad(Yd , params[bt].Mdup );           // (1/2) Yd^T Mdup Yd
double Emelt = 0.5 * quad(Ydm, params[bt].Mmelt) + params[bt].dG;

// --- soft-min mixing ---
double rho = params[bt].rho;
double m   = std::min(rho*Edup, rho*Emelt);                 // log-sum-exp stabilization
double ed  = std::exp(-(rho*Edup  - m));
double em  = std::exp(-(rho*Emelt - m));
double Zinv= 1.0/(ed+em);
double pdup = ed*Zinv, pmelt = em*Zinv;
double Else = -(1.0/rho)*(std::log(ed+em)) + (1.0/rho)*m;   // == -(1/rho) ln Z

// --- gradients dE/dYd ---
double gdup[6], gmelt[6], g[6];
matvec6(params[bt].Mdup , Yd , gdup );
matvec6(params[bt].Mmelt, Ydm, gmelt);
for (k=0;k<6;k++) g[k] = pdup*gdup[k] + pmelt*gmelt[k];

// --- barrier bump (twist only) ---
double db = phi_unw - params[bt].phi_b;
double bexp = params[bt].B * std::exp(-0.5*db*db/(params[bt].wb*params[bt].wb));
double Ebump = bexp;
g[2] += -db/(params[bt].wb*params[bt].wb) * bexp;          // dE_bump/dphi

// --- split into rotational A and translational B, then identical Y-mapping ---
double A[3] = { g[0], g[1], g[2] };
double Bv[3]= { g[3], g[4], g[5] };
// torque_2 = T1 Smat (Jinvtp A); force_1 = T1 Smat Bv; torque_1 = T1 (w x Smat Bv) + torque_2
// ... (verbatim from bond_rbp_fene.cpp lines 240-252) ...

if (eflag) bond_energy = Else + Ebump;   // replaces 0.5*Yd^T M Yd in ev_tally
```

`quad`, `matvec6` are trivial 6-D helpers. **Off-diagonal step couplings** (3-pt/4-pt) are handled
exactly as today: the angle/dihedral styles compute `Yd1^T M Yd2`. For the bimodal model these
coupling blocks are evaluated with the *softmax-selected* state (§6); to first order the recommended
implementation keeps the cross-couplings linear-harmonic (duplex) and adds cooperativity through a
separate junction term, see §6.

---

## 4. The tunable barrier (kinetics)

Three independent knobs, in increasing order of recommendation:

1. **`rho` (soft-min sharpness).** Raising `rho` sharpens the cusp where `E_dup = E_melt`, raising the
   intrinsic crossover ridge. But it couples to well geometry and can introduce stiff forces near the
   ridge. Keep `rho = 1` (physical Boltzmann sum) unless a sharper transition is explicitly wanted.

2. **`M_melt` magnitude (Sicard–Manghi mechanism).** Because the melted basin must be *traversed in
   twist* from `0` to `~-2pi`, a stiffer `M_melt` twist component raises the elastic cost of the
   intermediate (partially untwisted) configurations — the barrier scales affinely with
   `beta * kappa_phi^melt` (Sicard–Manghi 2015). This is the *physically motivated* barrier control:
   it ties the kinetic barrier to the same parameter that sets melted-state elasticity.

3. **Explicit Gaussian bump `B exp(-(phi-phi_b)^2/2w_b^2)` (RECOMMENDED for independent tuning).**
   Decouples kinetics from thermodynamics entirely: `B` sets barrier height, `phi_b` its location
   along twist, `w_b` the transition-state curvature, with **no effect on either basin minimum**
   (choose `phi_b` between basins and `w_b` small enough that the bump is negligible at both minima).
   Calibrate `B` so the total forward barrier `≈ 22 kBT` (10-bp bubble nucleation) and reverse
   `≈ 13 kBT` (closure) per Sicard–Manghi 2015 / Manghi–Destainville 2016; for buckling-dominated
   regimes target `~10 kBT` (Dittmore 2018). The asymmetry (forward vs reverse) is set automatically
   by `dG` plus the applied-torque tilt `-tau*phi`.

The barrier acts on the **unwrapped twist** and is therefore periodicity-free; it is the natural
place to encode the rate-limiting collective untwisting (Dasanna 2013: twist re-winding is the slow
coordinate).

---

## 5. Melted-state elasticity and the twist shift

### 5.1 Softening the stiffness (`M_melt`)

Build `M_melt` from the L-DNA / ssDNA elastic constants, coarse-grained to 10 bp with the same
Becker–Everaers / PolyCG dictionary used for `M_dup`:

- **Bending** (`M_melt[0][0]`, `M_melt[1][1]` — tilt, roll): persistence `A_melt ≈ 3-7 nm` (L-DNA,
  Sheinin–Wang 2011) down to `1-3 nm` (ssDNA, Maffeo 2014). At 10 bp/bead (`l = 3.4 nm` ds, `4.6 nm`
  L-DNA), the per-step bending stiffness `kappa_bend = A_melt / l_bead`, i.e. roughly an order of
  magnitude softer than duplex.
- **Twisting** (`M_melt[2][2]`): torsional modulus `C_melt ≈ 19-20 nm` (L-DNA) or `≈ 0` (fully
  denatured ssDNA). Soft twist is what lets the basin sit at `-2pi` cheaply.
- **Stretch** (`M_melt[5][5]`): softer; ssDNA is highly extensible. Keep the FENE non-extensibility
  cap (it already lives in `bond_rbpfene`) but lower the harmonic stretch.
- **Cross terms** (twist-stretch, twist-bend): set small; L-DNA twist-bend is poorly constrained, use
  `≈ 0` unless data dictates otherwise.

Concretely, a reasonable starting `M_melt = diag(A_melt/l, A_melt/l, C_melt/l, k_shift, k_shift, k_str)`
with `A_melt/l ≈ 0.1 * (A_dup/l)`, `C_melt/l ≈ 0.2 * (C_dup/l)`.

### 5.2 Shifting intrinsic twist by ~ -1 turn

`Delta[2] = dphi_melt`. Two ways to fix the number:

- **Exact -1 turn:** `dphi_melt = -2pi` per 10-bp bead. This is the "release one full negative turn
  on melting" picture used at 10-bp resolution.
- **L-DNA intrinsic twist:** L-DNA is ≈ -13 bp/turn = `-2pi/13` rad/bp. Over a 10-bp bead the *melted
  intrinsic* total twist is `10 * (-2pi/13) = -1.54 * 2pi`, vs duplex `10 * (2pi/10.5) = +0.95*2pi`.
  The basin-to-basin twist shift is then `dphi_melt = (-1.54 - 0.95)*2pi ≈ -2.49 * 2pi ≈ -15.6 rad`.
  The "~ -1 turn" headline is the *change in linking* (`ΔLk ≈ -1` per bead from the duplex value);
  the *intrinsic-twist* shift is larger because L-DNA is left-handed. **Make `dphi_melt` a per-bead
  database parameter** so either convention can be set from data; default to the L-DNA value.

### 5.3 Sequence-dependent offset `dG(seq)` from NN melting parameters

At 10-bp resolution each bead spans ~9 dinucleotide steps. Assign

```
dG_bead(seq) = dG_init + sum_{steps in bead} dG_NN(step) + dG_termAT + dG_salt          (Eq. 10)
```

using SantaLucia 1998 / SantaLucia–Hicks 2004 dG37 per step (kcal/mol, 1 M Na+):
GC -2.24, CG -2.17, GG -1.84, CA -1.45, GT -1.44, GA -1.30, CT -1.28, AA -1.00, AT -0.88, TA -0.58;
init +0.98/+1.03, symmetry +0.43. Convert with `1 kcal/mol ≈ 1.62 kBT` (300 K). Use
`dG(T) = dH - T dS` (full dH/dS tables) at the simulation temperature, then salt-correct
(Owczarzy 2004 monovalent; 2008 Mg2+). **Per-bead range ≈ +7 kBT (AT-rich) to +24 kBT (GC-rich)**,
i.e. duplex-favoring at zero torque. Cross-checks: Cocco–Marko per-bp `dG_AT≈1.1, dG_GC≈3.5 kBT`;
Marko–Neukirch sequence-averaged `eps_M ≈ 2.5 kBT/bp` (so ~25 kBT/bead) — note this includes the
elastic creation cost, part of which our model carries in `(1/2)Delta^T M_melt Delta`; therefore set
`dG` to the *configurational/H-bond* part and let the elastic shift supply the rest, then re-tune to
the Maxwell target (§5.4).

This `dG(seq)` is a **per-bond-type scalar** (the model already supports one bond type per position
for sequence-dependent fields), assigned by the database generator (`cgRBPTools`) from the bead's
10-bp window.

### 5.4 Calibration against the coexistence torque (the hard constraint)

The Maxwell / common-tangent construction (Marko 2007; Marko–Neukirch 2013) fixes the *relationship*
between offset and twist shift. Coexistence (degenerate basins under torque) requires

```
tau_c = dE_total / |dphi_melt|   where  dE_total = dG + (1/2)Delta^T M_melt Delta - (1/2)*0   (Eq. 11)
```

Set `tau_c = -10..-11 pN·nm ≈ -2.4..-2.7 kBT/rad` (since `kBT ≈ 4.1 pN·nm`, `tau[kBT/rad] =
tau[pN·nm]/4.1`). This is a **single scalar constraint** that pins the *sequence-averaged* `dG`
given `M_melt` and `dphi_melt`; sequence variation around the mean is then supplied by Eq. 10. Do not
import the oxDNA ~3 pN·nm plateau (model-specific). This calibration is the acceptance test in §8.

---

## 6. Cooperativity: melting domains / bubbles

Isolated per-step double wells give *no* cooperativity — every bead would melt independently. Real
melting is domain-like (bubbles) because creating a duplex/melted **junction** costs energy. Two
compatible routes, both reusing the existing banded `M`:

### 6.1 Domain-wall cost via the existing 3-pt/4-pt couplings (recommended, MD-native)

The off-diagonal coupling blocks already computed by `angle_rbp` / `dihedral_rbp`
(`E = Yd_i^T M_coup Yd_{i+1}`) are the natural carrier of a domain-wall penalty. Make the coupling
**state-dependent**: define a per-step "meltedness"

```
chi_i = p_melt(step i)      in [0,1]   (the softmax weight, already computed)
```

and add a junction term that penalizes *gradients* of `chi`:

```
beta*E_junction = sum_i  J * (chi_i - chi_{i+1})^2                                       (Eq. 12)
```

`J` is the domain-wall stiffness (Benham junction `a ≈ 10-11 kcal/mol ≈ 17 kBT`; tune to the ~22 kBT
nucleation target jointly with the bump `B`). Because `chi_i` is a smooth function of `Yd_i` (twist),
`dE_junction/dYd_i` is analytic:

```
dE_junction/dphi_i = 2 J [ (chi_i - chi_{i-1}) + (chi_i - chi_{i+1}) ] * dchi_i/dphi_i
dchi_i/dphi_i      = rho * chi_i (1 - chi_i) * (dEdup/dphi - dEmelt/dphi)|_i
```

This is a *next-nearest-neighbor* term, structurally identical to the existing 3-pt coupling, so it
slots into `angle_rbp` (which already reads `fwd_euler` for both sub-steps via `fix_rbp_lrf`). It
makes a sharp duplex/melted interface cost `~J`, producing finite bubbles and the constant-torque
plateau, and supplies the *nucleation* barrier complementary to the per-step `B`.

### 6.2 Transfer-matrix Ising layer (analysis / optional MC)

For equilibrium occupancy and the torque plateau (validation), reduce the model to the Storm–Nelson /
Benham 2-state transfer matrix `[[exp(-g_dup), exp(-J)],[exp(-J), exp(-(g_melt+dG-tau*dphi))]]` and
solve exactly (Fye–Benham 1999 handles the global fixed-Lk constraint). This is *not* run inside MD
but used to (a) set `J`, `dG`, `B` self-consistently and (b) provide the analytic reference curve the
MD must reproduce (§8). Optionally, a hybrid MC move flipping `chi_i -> {0,1}` can accelerate barrier
crossing, but is not required.

---

## 7. Implementation plan

### 7.1 New bond style `bond_rbpmelt` (mirror the FENE template)

Create `bond_rbp_melt.{h,cpp}` by cloning `bond_rbp_fene.{h,cpp}`. Changes to `RBPParams`:

```cpp
struct RBPParams {
   // --- existing duplex well (rename Mmat -> Mdup conceptually; keep blocks) ---
   double Ystatic[6]; double Smat[3][3]; double srot[3]; double svec[3];
   double Mdup[6][6];  double Mr[3][3], Mt[3][3], Mtr_bl[3][3], Mtr_tr[3][3];   // duplex blocks
   double equidist; bool subtract_groundstate;
   // --- FENE (kept; melted basin still needs the non-extensibility cap) ---
   double K, Rc, R0; bool fene_active; double Rspan, Rspan2;
   // --- NEW: melted basin ---
   double Mmelt[6][6];        // soft stiffness
   double Delta[6];           // ground-state shift; Delta[2]=dphi_melt, Delta[5]=drise_melt
   double dG;                 // sequence-dependent offset (kBT)
   double rho;                // soft-min sharpness (default 1)
   // --- NEW: barrier bump ---
   double Bbump, phi_b, w_b;
   // --- NEW: cooperativity ---
   double Jwall;              // domain-wall stiffness (used by angle/dihedral style)
   bool   melt_active;        // false => behaves exactly like bond_rbpfene
};
```

`compute()` gets the kernel of §3 (replacing only the wrench block; FENE and force application are
verbatim). When `melt_active == false` or `dG` huge, the soft-min collapses to the duplex well
(Eq. 7) — exact backward compatibility.

### 7.2 What must come from `fix_rbp_lrf` (Challenge 1)

Add to `FixRBPLRF` a persistent per-atom array `double *twist_unwrapped` (size `nmax`), keyed to the
left atom `id1` exactly like `fwd_euler`. In `compute_lrf` Phase 2, after computing the wrapped twist
`euler[2]`, update the winding counter:

```cpp
double prev   = twist_unwrapped[id1];          // last step's continuous twist (init from euler[2])
double wrapped= euler[2];                       // in (-pi, pi]
double base   = prev - std::remainder(prev, 2*M_PI);   // integer-turn part of prev
double cand   = base + wrapped;
// pick the 2pi branch closest to prev (unwrap across the frame)
if (cand - prev >  M_PI) cand -= 2*M_PI;
if (cand - prev < -M_PI) cand += 2*M_PI;
twist_unwrapped[id1] = cand;
```

This must be added to `grow_arrays`, `copy_arrays`, `set_arrays`, `pack_exchange`/`unpack_exchange`,
`pack_forward_comm`/`unpack_forward_comm` (extend the comm payload from 12 to 13 doubles), and — new —
**`write_restart`/`read_restart`** so the winding number survives restarts (the fix currently has
none; this is required for the unwrapped twist to be a true trajectory state variable). On first
step / `set_arrays`, initialize `twist_unwrapped[id1] = wrapped` (zero winding). Caveat (from the
review): per-frame unwrapping is unambiguous only if the per-step twist changes by `< pi` between
`pre_force` calls — guaranteed at 10-bp resolution with a reasonable timestep; assert/warn otherwise.

Cross-check (optional, analysis only): compute `Lk = Tw + Wr` per frame (Klenin–Langowski writhe over
bead centers) and verify `Tw = sum_i twist_unwrapped[i] / 2pi` is consistent — a topology audit, not
needed for forces.

### 7.3 Database / `parse_rbp.h` format changes

Add per-bond-type melted-basin coefficients after the existing `Ystatic(6) + M(48 or 21/12/6)` block.
Two clean options:

- **New bond style tag** `bond style: rbpmelt` (validated in `coeff` like the FENE check at
  `bond_rbp_fene.cpp:446`).
- **Per-type coefficient extension:** append `[ dG, dphi_melt, drise_melt, rho, Bbump, phi_b, w_b,
  Jwall, Mmelt(21 upper-tri) ]` = 28 extra coeffs. Add metadata flags `melt_set` (bool) and
  `melt_seqdep` (bool) mirroring `seqs_set`. The `cgRBPTools` database generator computes `dG(seq)`
  per bead (Eq. 10) and a (sequence-independent or sequence-scaled) `Mmelt`, `dphi_melt`.

Keep `subtract_groundstate`, `seqs_set`, `closed`, `unit_length` semantics unchanged. `Mmelt` is read
and block-decomposed by the same `assign_blocks_` machinery (add `Mmelt_r/_t/_tr_*` blocks).

### 7.4 Applied-torque / fixed-Lk ensembles

- **Constant torque (Gibbs):** add `-tau * phi_unw(i)` to each step energy; its gradient is a constant
  torque `-tau` about the twist axis, applied through the same `Jinvtp` mapping (a one-line addition
  to `g[2]`). Or, more cleanly, apply `fix addtorque`/`fix addtorque/atom` to the terminal beads
  (LAMMPS native) and clamp the other end — no force-field change. The melting torque makes the melted
  basin the global minimum.
- **Fixed Lk:** clamp both ends (`fix rigid` or harmonic angular springs) and inject `n` turns at
  setup; `Lk = Tw + Wr` is then conserved and `tau` is *measured* as the conjugate response
  (Bryant RBT-style). The unwrapped twist (§7.2) + Klenin–Langowski writhe give `Tw` and `Wr`.

### 7.5 Collective variables for enhanced sampling

Expose `phi_unw` (sum over a region) and the melted-fraction `sum_i chi_i` as per-atom/global compute
outputs so PLUMED (`fix plumed`) can bias them. **Do not use PLUMED's built-in TORSION CV** (wrapped
to `(-pi,pi]`); feed the unwrapped accumulator via a `CUSTOM`/`MATHEVAL` CV, or bias the
periodicity-free melted-step count. Then: metadynamics / umbrella+WHAM/MBAR for the 2-D PMF
`(Tw_unw, melted-fraction)`, FFS for rates with the melted-step count as order parameter, Kramers /
string method for the 1-D reduced quartic.

---

## 8. Validation plan

1. **Backward-compatibility (must pass first).** Set `melt_active = false` (or `dG = +1e6`). The new
   style must reproduce `bond_rbpfene` energies and forces to machine precision on an existing test
   system (duplex elastic moduli A, C from equilibrium fluctuation covariance must match the input
   `Mdup`). Confirms Eq. 7 numerically.

2. **Single-step double-well shape.** One step, scan imposed twist `phi` from `+0.5` to `-3` turns
   (umbrella on `phi_unw`), reconstruct `G(phi)` by WHAM. Verify: two minima at `0` and `dphi_melt`,
   offset `dG`, barrier height `= B (+ junction)`, and that curvatures match `Mdup[2][2]`,
   `Mmelt[2][2]`. Cross-check against the analytic quartic fit (Eq. 5).

3. **Coexistence torque (the key physics target).** Apply constant torque `tau` to a homogeneous
   chain; sweep `tau`; measure equilibrium melted fraction `<chi>`. The transition must be a
   **constant-torque plateau** at `tau_c = -10..-11 pN·nm` (Marko 2007; Sheinin–Wang 2011). Tune the
   sequence-averaged `dG` (Eq. 11) until `tau_c` matches. Two-state occupancy vs torque must follow
   the transfer-matrix/Ising prediction (§6.2).

4. **Extension/twist vs turns ("hat curve").** Fixed-Lk ensemble: clamp ends, inject `n` turns, plot
   extension and torque vs `n`. Below `F_c≈0.5 pN` torque rises linearly then plateaus; above it,
   negative turns are absorbed by melted bubbles at constant extension (Strick 1998; Allemand 1998).

5. **Sequence dependence.** AT-rich vs GC-rich beads: bubbles must nucleate preferentially in AT-rich
   regions (Matek 2015: 84-91% AT in bubble tips), and `tau_c` / `Tm` ordering must follow the
   SantaLucia `dG(seq)` ordering.

6. **Cooperativity / bubble size.** Vary `Jwall`: larger `J` => fewer, larger bubbles (sharper
   transition); `J=0` => independent per-step melting (broad). Compare bubble-size distribution to
   Poland–Scheraga (cooperativity `sigma_coop ~ 1e-4..1e-5`).

7. **Kramers rates / hysteresis.** From the `G(phi)` barrier + a measured diffusion coefficient,
   predict opening/closing rates (Kramers); compare to direct FFS rates and to torque-jump hysteresis
   loops (finite barrier => hysteresis; Lebel 2014 dwell-time kinetics; target opening barrier
   ~22 kBT, closure ~13 kBT, Sicard–Manghi 2015).

8. **Topology audit.** Throughout, verify `Lk = Tw + Wr` holds each frame (Klenin–Langowski writhe +
   unwrapped twist sum), confirming the released turns are correctly partitioned between local
   untwist (bubbles) and global writhe.

---

## 9. Summary of recommendation

- **Energy:** per-step **log-sum-exp of two shifted Gaussian wells** (Eq. 3/6) with `rho = 1`
  (KWLC/Wiggins Boltzmann sum), a stiff duplex well (existing `M`, origin in Y convention) and a
  **soft melted well** (`M_melt` from L-DNA/ssDNA, §5.1) shifted by `dphi_melt ≈ -2pi` in the
  **unwrapped** twist and raised by sequence-dependent `dG(seq)` (SantaLucia NN, Eq. 10).
- **Barrier:** an explicit Gaussian **bump on the unwrapped twist** (Eq. 9) for independent kinetic
  tuning, physically backed by the `M_melt` torsional-stiffness mechanism (Sicard–Manghi).
- **Cooperativity:** a **domain-wall penalty `J(chi_i - chi_{i+1})^2`** carried by the existing
  3-pt/4-pt couplings (Eq. 12), giving finite bubbles and the constant-torque plateau.
- **Exact reduction:** with the melted basin disabled it is **bit-for-bit the current harmonic model**
  (Eq. 7) — strict superset, safe to merge.
- **Gradient/forces:** softmax-weighted sum of the two harmonic gradients (Eq. 8) plus the bump
  (Eq. 9), mapped to body torques by the **existing `leftJacobianInverseTransposed` / Smat / T1**
  pipeline — only the (A, B) vectors change.
- **Implementation:** new `bond_rbpmelt` style cloned from `bond_rbp_fene`; `fix_rbp_lrf` extended with
  a **persistent unwrapped-twist winding counter + restart** (Challenge 1) and a 13th comm double;
  database/`parse_rbp.h` extended with melted-basin coefficients and `melt_set`/`melt_seqdep` flags.
- **Calibration:** sequence-averaged `dG` pinned by the Maxwell construction to `tau_c ≈ -10 pN·nm`;
  `B`+`J` tuned to ~22 kBT nucleation / ~13 kBT closure.
- **Validation:** backward-compat → single-step PMF → constant-torque plateau → hat curve → sequence
  dependence → cooperativity → Kramers rates, with a per-frame `Lk = Tw + Wr` audit.
