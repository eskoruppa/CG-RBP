# CG-RBP — Coarse-grained Rigid Base Pair force field for LAMMPS

A LAMMPS user package implementing a **sequence-dependent, variable-resolution
coarse-grained rigid base pair (RBP) model of double-stranded DNA**. Each DNA
base pair (or a coarse-grained block of several base pairs) is represented by an
oriented rigid body — a bead carrying a position **and** an orthonormal triad —
and the sequence-dependent elastic energy is realised through custom bond,
angle, and dihedral styles that act on the relative pose of neighbouring beads.

The package provides the C++ interaction styles only. Model setup (sequence →
topology + parameters), trajectory analysis, back-mapping, and visualisation are
handled by the companion Python package
[**cgRBPTools**](https://github.com/eskoruppa/cgRBPTools) (see
[Companion tooling](#companion-tooling-cgrbptools)).

> Contributing author: **Enrico Skoruppa** (School of Physics and Astronomy,
> University of Edinburgh). Distributed under the GNU GPL, consistent with
> LAMMPS. See [Citation](#citation).

---

## Contents

- [The model in brief](#the-model-in-brief)
- [Package contents](#package-contents)
- [Installation](#installation)
- [Interaction styles](#interaction-styles)
- [Simulation requirements](#simulation-requirements)
- [Quick start](#quick-start)
- [Specifying coefficients](#specifying-coefficients)
- [Ground-state convention (X vs Y)](#ground-state-convention-x-vs-y)
- [The auxiliary fix `rbp/lrf`](#the-auxiliary-fix-rbplrf)
- [Data file format](#data-file-format)
- [Database (`.db`) file format](#database-db-file-format)
- [The FENE bond (`rbpfene`)](#the-fene-bond-rbpfene)
- [Output and analysis](#output-and-analysis)
- [Companion tooling: cgRBPTools](#companion-tooling-cgrbptools)
- [Units, restarts, and other notes](#units-restarts-and-other-notes)
- [Citation](#citation)

---

## The model in brief

Each base pair is a rigid frame $\tau_i \in \mathrm{SE}(3)$ — a triad
$\mathcal{T}_i \in \mathrm{SO}(3)$ (orientation) together with a position
$\mathbf{r}_i \in \mathbb{R}^3$. The junction between adjacent beads is the
body-frame step transformation $g_i \equiv \tau_i^{-1}\tau_{i+1} \in \mathrm{SE}(3)$,
whose six coordinates $\mathbf{X}_i = (\mathbf{\Omega}_i,\ \mathbf{w}_i)^\intercal \equiv \mathcal{P}(g_i)$
comprise three rotational components $\mathbf{\Omega}_i$ (tilt, roll, twist) and
three translational components $\mathbf{w}_i$ (shift, slide, rise) —
**rotational components first**.

Sequence dependence enters through a ground-state shape $\mathbf{X}_{0i}$ and a
Gaussian penalty on the deformation $\mathbf{X}_{\Delta i} = \mathbf{X}_i - \mathbf{X}_{0i}$.
Collecting all step deformations into $\bar{\mathbf{X}}_\Delta$, the elastic
energy is

```math
\beta E = \tfrac{1}{2}\,\bar{\mathbf{X}}_\Delta^{\intercal}\, M\, \bar{\mathbf{X}}_\Delta ,
```

with $\beta = 1/k_\mathrm{B}T$ and $M$ a banded stiffness matrix. The package
realises the 6×6 blocks of $M$ as standard LAMMPS bonded interactions, according
to how many beads each coupling involves. The diagonal, first-, and
higher-off-diagonal blocks contribute

```math
\begin{aligned}
\beta E^{(0)}_i &= \tfrac{1}{2}\,\mathbf{X}_{\Delta i}^{\intercal} M^{(0)}_i \mathbf{X}_{\Delta i}, \\
\beta E^{(1)}_i &= \mathbf{X}_{\Delta i}^{\intercal} M^{(1)}_i \mathbf{X}_{\Delta,i+1}, \\
\beta E^{(m)}_i &= \mathbf{X}_{\Delta i}^{\intercal} M^{(m)}_i \mathbf{X}_{\Delta,i+m}\quad (m \ge 2),
\end{aligned}
```

and map onto bonds, angles, and dihedrals respectively:

| Block | Coupling | LAMMPS object | Style(s) |
|---|---|---|---|
| M⁽⁰⁾ (diagonal) | local, 2 beads | **bond** | `rbp`, `rbpfene` |
| M⁽¹⁾ (first off-diagonal) | range-1, 3 beads | **angle** | `rbp` |
| M⁽ᵐ⁾, m ≥ 2 (higher) | range-m, 4 beads | **dihedral** | `rbp` |

A local coupling relates the two triads of one junction → a **bond**. A
range-$m$ coupling correlates two junctions, $\mathbf{X}_{\Delta i}$ and
$\mathbf{X}_{\Delta,i+m}$; for $m = 1$ they share the central triad (three beads
→ an **angle**), and for $m \ge 2$ all four beads are distinct (a **dihedral**
over beads $i, i+1, i+m, i+m+1$). Because an angle or dihedral carries two
junctions, its coefficients specify **two** ground-state vectors and the **full
(non-symmetric) 6×6 coupling block**.

The elastic torques on the oriented beads are evaluated in the Lie-algebra
"force-wrench" formulation; the per-bead triad and inverse-transposed left
Jacobian needed for this are precomputed once per timestep by the auxiliary
[`fix rbp/lrf`](#the-auxiliary-fix-rbplrf).

---

## Package contents

| File | Purpose |
|------|---------|
| `bond_rbp.cpp/.h`      | Bond style `rbp` — local (diagonal) block M⁽⁰⁾ |
| `bond_rbp_fene.cpp/.h` | Bond style `rbpfene` — `rbp` plus a built-in one-sided FENE non-extensibility term |
| `angle_rbp.cpp/.h`     | Angle style `rbp` — range-1 coupling block M⁽¹⁾ |
| `dihedral_rbp.cpp/.h`  | Dihedral style `rbp` — range-m coupling blocks M⁽ᵐ⁾, m ≥ 2 |
| `fix_rbp_lrf.cpp/.h`   | Auxiliary fix `rbp/lrf` — per-timestep triads / local reference frames (auto-created) |
| `parse_rbp.cpp/.h`     | Reader for the `.db` parameter database (metadata + coefficient sections) |
| `so3.h`                | Header-only SO(3)/SE(3) utilities (Rodrigues exp/log, left Jacobian, triad↔euler) |
| `lamath.h`             | Header-only small 3×3 / 3-vector linear algebra and positive-definiteness (Cholesky) checks |
| `install.sh`           | LAMMPS package install/uninstall hook (checks MOLECULE + ASPHERE) |

The companion Python package, the manuscript sources, and research notes
(`cgRBPTools/`, `Manuscript_cgRBP/`, `melting_design/`, `oxDNA/`) are **not part
of the tracked package** — they are excluded via `.gitignore`.

---

## Installation

CG-RBP is an optional LAMMPS package. It **requires** the `MOLECULE` package
(bond/angle/dihedral topology) and the `ASPHERE` package (ellipsoidal particles
and rigid-body integrators). Both are mandatory and are checked by `install.sh`.

From the LAMMPS `src/` directory:

```bash
# copy CG-RBP into src/ (if not already present)
# then enable it together with its dependencies:
make yes-molecule yes-asphere yes-cg-rbp

# build as usual, e.g. with MPI:
make mpi
```

The CMake build works the same way — enable `PKG_MOLECULE`, `PKG_ASPHERE`, and
`PKG_CG-RBP`.

---

## Interaction styles

The package registers the following styles:

- `bond_style rbp`        — local block, harmonic in the step deformation
- `bond_style rbpfene`    — `rbp` + one-sided shifted FENE (non-extensibility)
- `angle_style rbp`       — range-1 (next-neighbour junction) coupling
- `dihedral_style rbp`    — a single style representing **all** range-$m \ge 2$
  couplings
- `fix rbp/lrf`           — auxiliary precompute fix, **created automatically**;
  users neither add nor configure it (see below)

All bonded styles operate on oriented rigid bodies and therefore enforce, at
initialisation:

- an `atom_style` carrying an orientation (ellipsoids) — see below;
- a three-dimensional simulation;
- bead orientations evolving under a **rotational** integrator
  (`fix nve/asphere`); a plain `fix nve` leaves the quaternions frozen.

---

## Simulation requirements

| Requirement | Reason |
|-------------|--------|
| `atom_style hybrid molecular ellipsoid` | beads need topology **and** an orientation (quaternion) |
| `dimension 3` | the styles abort in 2D |
| `fix nve/asphere` (+ `fix langevin ... angmom`) | integrates both positions and quaternions; thermostats translational and rotational DOF |
| `units lj` (recommended) | parameters are in units of `k_BT` and a chosen length unit |
| consistent 5′→3′ bond orientation | step quantities are keyed to the **left** bead of each step (see [fix `rbp/lrf`](#the-auxiliary-fix-rbplrf)) |

No `pair_style` is required for a bare elastic chain. Excluded-volume or
electrostatic interactions may be added with any standard pair style. If you do
add a pair style, keep the **1–2 special-bonds weight at zero**
(`special_bonds lj/coul 0 x x`); the styles emit a warning if `special_lj[1] ≠ 0`.

---

## Quick start

### Minimal homopolymer (inline coefficients, local couplings only)

A chain whose stiffness matrix is purely local needs a single bond type and no
angle/dihedral terms:

```lammps
units           lj
dimension       3
boundary        p p p
atom_style      hybrid molecular ellipsoid

read_data       model.data

bond_style      rbpfene
#          type  K    Rc   R0    |  ground state Y0 (rot,trans)  |  6 diagonal moduli
bond_coeff  *   500  1.15 1.3  &
            0.0 0.0 0.0 0.0 0.0 1.0 &
            11.7647 11.7647 29.41 200 200 200

fix     1 all nve/asphere
fix     2 all langevin 1.0 1.0 1.0 12345 angmom 3.0

timestep        0.005
thermo          1000
run             1000000
```

Here $\mathbf{Y}_0 = (0,0,0,0,0,1)^\intercal$ is a straight, untwisted step of
unit rise; the six diagonal moduli are two bending, one twist
(`11.76, 11.76, 29.41`) and three translational (`200`).

### Sequence-dependent chain (parameters from a database file)

For a real sequence, every position has its own parameter set, so the bonded
*types* enumerate positions along the chain. Parameters are read from a `.db`
database and the matching `.data` file — both generated by
[cgRBPTools](#companion-tooling-cgrbptools). Example for an open chain of 271
beads with couplings up to range `m = 2`:

```lammps
units           lj
dimension       3
boundary        p p p
atom_style      hybrid molecular ellipsoid

read_data       model.data

# elastic interactions read from the database file
bond_style      rbpfene
bond_coeff      1*270 dbfile model.db 1
angle_style     rbp
angle_coeff     1*269 dbfile model.db 1
dihedral_style  rbp
dihedral_coeff  1*268 dbfile model.db 1

# translational + rotational integration
fix     1 all nve/asphere
fix     2 all langevin 1.0 1.0 1.0 12345 angmom 3.0

timestep        0.005
thermo          1000
run             1000000
```

For an open chain of $N+1$ beads with maximal range $m_\mathrm{max}$, there are
$N$ bond types, $N-1$ angle types, and $\sum_{m=2}^{m_\mathrm{max}}(N-m)$
dihedral types (for the 271-bead, $m_\mathrm{max}=2$ example: 270 / 269 / 268).
Closed molecules wrap the couplings around the seam and the counts increase
accordingly.

---

## Specifying coefficients

Each type may be specified **inline** or read **from a database file**. In every
case the six-vectors list **rotational components before translational**, and
stiffness entries are given in **row-major** order.

### Inline

| Style | Argument count | Layout |
|-------|----------------|--------|
| `bond rbp`     | **12 / 18 / 27** | 6 ground state + M⁽⁰⁾ as 6 diagonal, 12 block-diagonal, or 21 upper-triangular entries |
| `bond rbpfene` | **15 / 21 / 30** | 3 FENE (`K Rc R0`) + the 12/18/27 above |
| `angle rbp`    | **48** | 6 + 6 ground states + full 6×6 block (36 entries) |
| `dihedral rbp` | **48** | 6 + 6 ground states + full 6×6 block (36 entries) |

The three `M` layouts for bonds are: **diagonal** (6), **block-diagonal**
(upper-triangular within each 3×3 rotation/translation block, 12), or **full**
symmetric (upper-triangular of the whole 6×6, 21). The full local block is
checked to be **positive definite**; a non-PD block aborts the run. Angle and
dihedral coupling blocks are in general **not symmetric**, so all 36 entries are
required and no PD check is applied.

### From a database file

```lammps
<bond|angle|dihedral>_coeff <type-range> dbfile <file> <start-id>
```

The database index runs over the type range starting at `<start-id>` and
incrementing by one; e.g. `bond_coeff 1*270 dbfile model.db 1` maps bond type
`i` to database bond entry `i`. The style named in the script is checked against
the corresponding `... style` field in the database metadata, and the convention
(`subtract groundstate`) is taken from the metadata.

---

## Ground-state convention (X vs Y)

How the sequence-dependent ground state is removed from the energy is controlled
by the boolean `subtract_groundstate` (metadata field `subtract groundstate`:
`1 = X`, `0 = Y`):

- **Y convention (default, `0`/`false`) — recommended.** The ground state is
  factored out at the group level, $g_i = s_i d_i$, and the energy
  $\beta E_\mathrm{Y} = \tfrac{1}{2}\,\bar{\mathbf{Y}}_\Delta^{\intercal} M_\mathrm{Y}\,\bar{\mathbf{Y}}_\Delta$
  is quadratic in the dynamic coordinate $\mathbf{Y}_{\Delta i}$. This is
  numerically robust for arbitrarily large intrinsic twist (including near
  180°).
- **X convention (`1`/`true`).** The ground state is subtracted additively in
  coordinate space, $\mathbf{X}_{\Delta i} = \mathbf{X}_i - \mathbf{X}_{0i}$,
  giving $\beta E_\mathrm{X} = \tfrac{1}{2}\,\bar{\mathbf{X}}_\Delta^{\intercal} M_\mathrm{X}\,\bar{\mathbf{X}}_\Delta$.

Inline coefficients always use the **Y** convention. **All interaction types
sharing a junction must use the same convention** — this is enforced at
initialisation (see below).

---

## The auxiliary fix `rbp/lrf`

The first time any `rbp`/`rbpfene` style initialises, the package automatically
creates `fix rbp_lrf all rbp/lrf`. It:

1. **Precomputes, once per timestep**, the triad of every bead (local + ghost)
   and, per bond, the step's euler vector and the inverse-transposed left
   Jacobian $J_L^{-\intercal}$ used to convert generalised forces to torques. The
   interaction styles read these shared quantities by pointer, avoiding
   redundant quaternion→matrix work.
2. **Validates junction consistency** at setup: all bond, angle, and dihedral
   types sharing a junction must declare the **same convention and the same
   ground-state** (`srot`, `svec`); a mismatch aborts the run.

Users neither add nor configure this fix, but it appears in the fix list and in
restart files.

**Important structural requirements** (enforced, but not obvious — relevant if
you build topologies by hand rather than with cgRBPTools):

- **Every angle/dihedral step must also carry an RBP bond.** The precompute fix
  sources its triad/euler/Jacobian data from the *bond* list; an angle or
  dihedral sitting on a step without a matching RBP bond has no data source and
  is rejected with a hard error.
- **Bonds, angles, and dihedrals must share a consistent head→tail (5′→3′)
  orientation.** Step quantities are defined relative to the **left** bead of a
  step and are *not* invariant under swapping the two beads, so orientation is
  physically meaningful. A reversed bond is caught as an error.
- Each bead must be the **left endpoint of at most one RBP bond** (holds for a
  linear backbone).

---

## Data file format

Read with `read_data`. Uses the `atom_style hybrid molecular ellipsoid` column
order — after `id type x y z` come the molecular sub-style's `mol` and the
ellipsoid sub-style's `ellipsoidflag density`:

```
4 atoms
3 bonds
2 angles
1 dihedrals
1 atom types
4 bond types
3 angle types
2 dihedral types
4 ellipsoids

-50.0 50.0 xlo xhi
-50.0 50.0 ylo yhi
-50.0 50.0 zlo zhi

Masses

1 1.0

Atoms      # id type x y z mol ellipsoidflag density

1 1 0.0 0.0 0.0 1 1 1
2 1 0.0 0.0 1.0 1 1 1
...

Ellipsoids # id shape_x shape_y shape_z qw qi qj qk

1 1.0 1.0 1.0 1.0 0.0 0.0 0.0
2 1.0 1.0 1.0 0.9553 0.0 0.0 0.2955
...

Bonds      # index type atom1 atom2
1 1 1 2
...

Angles     # index type atom1 atom2 atom3
1 1 1 2 3
...

Dihedrals  # index type atom1 atom2 atom3 atom4
1 1 1 2 3 4
```

Each bead carries `ellipsoidflag = 1` and a near-spherical shape (which sets the
moment of inertia for the rigid-body integrator) plus a unit quaternion defining
the bead triad. Connectivity lists two beads for a local junction (bond), three
for a range-1 coupling (angle), and the four beads `i, i+1, i+m, i+m+1` for a
range-`m` coupling (dihedral). The `type` numbers must match the database
identifiers assigned by the `*_coeff` commands.

---

## Database (`.db`) file format

The parameter database is a plain-text file read by the `dbfile` keyword. It is
normally produced by cgRBPTools as a **single combined file** that also carries
the sequence and connectivity used to build the `.data` file. LAMMPS reads
**only** the metadata block and the three coefficient sections (`Bond Coeffs`,
`Angle Coeffs`, `Dihedral Coeffs`); the `Seqs`, `Bonds`, `Angles`, and
`Dihedrals` sections are silently ignored on the LAMMPS side.

> **Note on the file extension.** The `dbfile` keyword accepts **any** filename;
> the extension is not checked. cgRBPTools writes `.db`.

### Metadata (leading `key: value` lines)

In the shipped configuration (`parse_rbp.h`) the following fields are **required**
— a missing field aborts the run. Coefficient identifiers must start at `1` and
be contiguous, and the counts/styles are cross-checked against the actual data.

| Field | Meaning |
|-------|---------|
| `number of rigid bodies` | number of beads in the molecule |
| `coupling range` | maximal range `m` (`0` = local only) |
| `number of bonds` / `angles` / `dihedrals` | number of terms of each kind |
| `number of bond types` / `angle types` / `dihedral types` | number of distinct parameter sets |
| `chars per atom` | base pairs per coarse-grained bead |
| `bond style` / `angle style` / `dihedral style` | expected style (`rbp` / `rbpfene`) |
| `subtract groundstate` | convention: `1 = X`, `0 = Y` |
| `seqs set` | whether per-bead sequence data is present (`1`/`0`) |
| `seqs centered` | whether retained triads are centred in each CG block (`1`/`0`) |
| `closed` | circular molecule (`1`/`0`) |
| `unit length` | length unit (nm) used during generation |

The only **optional** field is `unit energy` (energy unit in `k_BT`, default
`1.0`). Unknown metadata keys are ignored for forward compatibility.

> ⚠️ **Documentation note.** This required-field list is authoritative to the
> parser and is broader than the table in the manuscript SI: `seqs set`,
> `seqs centered`, `closed`, and `unit length` are **mandatory** in this build
> (only `unit energy` is optional). A database written to match the manuscript
> table alone would abort with e.g. *"Missing required metadata 'seqs set'"*.
> Databases produced by cgRBPTools include all of these fields and load cleanly.
> To relax a requirement, comment out the corresponding `RBP_META_REQUIRE_*`
> define in `parse_rbp.h` and recompile.

### Coefficient sections

Each line begins with an integer identifier and lists coefficients in the same
order as the inline form:

- **`Bond Coeffs`** — `id` + 6 ground state + stiffness (6 / 12 / 21 entries; a
  full row-major 6×6 = 36 is also accepted by the parser).
- **`Angle Coeffs` / `Dihedral Coeffs`** — `id` + 6 + 6 ground states + full
  6×6 block (36 entries, row-major).

---

## The FENE bond (`rbpfene`)

`bond_style rbpfene` adds, on top of the harmonic RBP wrench, a **one-sided,
shifted FENE** term that acts as a non-extensibility wall — silent for short
bonds and diverging as the bead separation approaches a limit. It is useful to
prevent chain crossing / overstretch without perturbing the equilibrium
elasticity.

The three FENE parameters `K Rc R0` precede the ground-state + stiffness entries
in the coefficient list. Writing $\Delta = R_0 - R_c$, $k_\mathrm{eff} = K\Delta^2$,
and the log argument $\rho(r) = 1 - (r-R_c)^2/\Delta^2$,

```math
E_\mathrm{FENE}(r) =
\begin{cases}
0, & r < R_c \\
-\tfrac{1}{2}\, k_\mathrm{eff}\, \ln \rho(r), & R_c \le r < r_\mathrm{max} \\
E(r_\mathrm{max}) + F(r_\mathrm{max})\,(r - r_\mathrm{max}), & r \ge r_\mathrm{max}
\end{cases}
```

- **$R_c$** — onset distance. Below $R_c$ the term is silent ($E = F = 0$), so
  the harmonic RBP elasticity governs small fluctuations; the potential is
  $C^1$-continuous at $R_c$.
- **$R_0$** — divergence radius; $E, F \to \infty$ as $r \to R_0$.
- **$K$** — bare stiffness. Setting $R_c = 0$ recovers the standard
  Kremer–Grest FENE.

**Robustness.** To survive integrator overshoot, the log argument is capped at
$\rho_\mathrm{min} = 0.1$ (i.e. $r_\mathrm{max} = R_c + \Delta\sqrt{1-\rho_\mathrm{min}}$);
beyond $r_\mathrm{max}$ the force is frozen and the energy linearly extrapolated.
A throttled warning is logged when $\rho < \rho_\mathrm{min}$, and a hard error is
raised if the bond is far past breaking ($\rho \le -3$). FENE is
**auto-deactivated** (with a warning) if $K \le 0$ or $R_0 \le R_c$; a passed
$R_c < 0$ is clamped to zero. When inactive, `rbpfene` behaves exactly like `rbp`.

See the manuscript appendix (`SISections/fene.tex`) for the full derivation.

---

## Output and analysis

Because each bead carries an orientation, a trajectory is fully specified only by
positions **and** quaternions. Expose the quaternions via a compute and write
them in a custom dump:

```lammps
compute q all property/atom quatw quati quatj quatk
dump    1 all custom 1000 traj.dump &
        id type x y z ix iy iz &
        c_q[1] c_q[2] c_q[3] c_q[4]
dump_modify 1 sort id
```

Include image flags (`ix iy iz`) so trajectories can be unwrapped in
post-processing. Elastic and kinetic energies come from standard computes; the
rigid-body rotational kinetic energy is reported separately:

```lammps
compute  erot all erotate/asphere
compute  epot all pe
thermo_style custom step temp pe c_erot etotal
```

Conversion of dumped quaternions into triads and positions, evaluation of mean
shapes and stiffnesses, and visualisation (including rendering bead triads in
ChimeraX) are handled by cgRBPTools.

---

## Companion tooling: cgRBPTools

Model setup and analysis live in the separate Python package
[**cgRBPTools**](https://github.com/eskoruppa/cgRBPTools) — *not* bundled in this
repository. It provides a complete sequence → simulation → analysis workflow:

```bash
git clone --recurse-submodules -j8 git@github.com:eskoruppa/cgRBPTools.git
cd cgRBPTools
pip install -e .          # requires Python ≥ 3.9; installs numpy, scipy, numba, matplotlib
```

The `--recurse-submodules` flag is required to fetch the **PolyCG** library and
the **cgNA+** parameter database.

**Generate a database + data file** for a DNA sequence, coarse-graining to one
bead per 10 bp, keeping couplings up to range `m = 2`, folding in a FENE term,
and writing a ground-state initial configuration:

```bash
python -m cgrbptools.lmp_input \
    -seqfn sequence.seq -m cgnaplus \
    -cg 10 -cr 2 -fene 200 1.1 1.35 \
    -conf gs -o model
```

Principal options (see the cgRBPTools README for the full list):

| Option | Meaning |
|--------|---------|
| `-seqfn` / `-seq` | DNA sequence (file or string) |
| `-m` | parameter model: `cgnaplus` (default), `md`, or `crystal` |
| `-cg` | base pairs per bead (default `1`) |
| `-cr` | coupling range `m_max` (default `1`; `0` = local only) |
| `-fene K Rc R0` | include a built-in FENE bond |
| `-closed` | generate a circular molecule (sequence length must be a multiple of `-cg`) |
| `-conf` | write an initial configuration: `straight`/`str`, `circular`/`circ`, `ground_state`/`gs` |
| `-dlk` | excess linking number for circular / supercoiled configurations |
| `-ul` / `-ue` | length (nm) / energy (`k_BT`) unit rescaling |
| `-vis` / `-pdb` / `-xyz` | also write ChimeraX / PDB / XYZ visualisations |
| `-o` | output base name |

Other modules: `cgrbptools.validate` (compare simulated mean shape and stiffness
against the reference parameters), `cgrbptools.io.parse_custom` (load custom
dumps → SE(3) poses), `cgrbptools.core.backmap` (reconstruct base-pair
resolution from coarse beads), and `cgrbptools.visualize`.

> Parameter models are supplied by the PolyCG backend: `cgnaplus` (Sharma et
> al. 2023, default), `md` (Lankaš et al. 2003), and `crystal` (Olson et
> al. 1998).

---

## Units, restarts, and other notes

- **Units.** Parameters are expressed in `k_BT` and a chosen length unit, so
  `units lj` is the natural choice. The `unit length` / `unit energy` metadata
  fields are recorded by the generator for bookkeeping; rescaling happens in
  cgRBPTools at generation time — LAMMPS does **not** re-scale coefficients by
  them.
- **Restarts.** All bonded styles are restart-capable: `write_restart` stores
  the primary state (ground state, `M`, convention; plus `K, Rc, R0` for FENE)
  and derived quantities are recomputed on read. The auxiliary `fix rbp/lrf`
  appears in restart files.
- **`write_data`.** The styles can emit their coefficients into a LAMMPS data
  file (`id` + 6 ground state + 21 upper-triangular `M` entries for bonds).
- **Tuning the parser.** Metadata requirements are compile-time switches
  (`RBP_META_REQUIRE_*` in `parse_rbp.h`); the precompute path can be disabled
  per style via the `*_PRECOMPUTE_ACTIVE` defines in the style headers (falling
  back to a slower, bond-independent local computation).

---

## Citation

If you use this package, please cite the accompanying manuscript,

> E. Skoruppa, D. Marenduzzo, and O. Henrich, *Molecular dynamics force fields
> for variable resolution sequence-dependent modelling of double-stranded DNA*
> (in preparation),

and the coarse-graining methodology it builds on:

> E. Skoruppa and H. Schiessel, *Systematic coarse-graining of sequence-dependent
> structure and elasticity of double-stranded DNA*, Phys. Rev. Research **7**,
> 013044 (2025). [doi:10.1103/PhysRevResearch.7.013044](https://doi.org/10.1103/PhysRevResearch.7.013044)

Depending on the parameter set used, please also cite the underlying elasticity
database: cgNA+ (Sharma et al. 2023), Lankaš et al. (2003), or Olson et al.
(1998).

## License

GNU General Public License, consistent with LAMMPS. See `LICENSE`.
