# Torque-Induced DNA Melting as Barrier Crossing in a Coarse-Grained Rigid-Base-Pair Model: A Master Literature Review

**Purpose.** This review consolidates 16 literature surveys into a single annotated reference supporting the extension of an existing SE(3) coarse-grained rigid-base-pair (CG-RBP) DNA force field (single harmonic well, "Y convention" ground-state factoring, banded 2-pt/3-pt/4-pt couplings) to model **torque-induced duplex melting at ~10 bp / bead resolution, treated as a barrier-crossing problem.** Two hard problems frame everything: **(1)** tracking twist beyond ±180° (a rotation matrix cannot store a winding number), and **(2)** building a phenomenological **bimodal** (stiff duplex / soft melted) free energy with a sequence-dependent ΔG offset and a tunable kinetic barrier.

A note on units used throughout: 1 kcal/mol ≈ 1.69 kBT at 37 °C ≈ 1.62 kBT at 300 K; kBT ≈ 4.1 pN·nm at 300 K. Persistence lengths are quoted in nm (a stiffness divided by kBT).

---

## 1. Experimental phenomenology of torque-induced melting

Key numbers established by this literature, which any model must reproduce:
- **Melting (denaturation) torque plateau: τ_melt ≈ −10 to −11 pN·nm** (negative = underwinding). This is the single most-cited target number, confirmed by rotor-bead tracking, angular optical traps, and magnetic torque tweezers.
- **Melted ("L-DNA") basin elasticity:** left-handed intrinsic twist ≈ −13 bp/turn; bending persistence ≈ 3 nm (vs ~43–50 nm for B-DNA); torsional modulus ≈ 20 nm (vs ~100 nm); rise 0.48 nm/bp (~40% longer than 0.34 nm/bp).
- **Supercoiling density at melting:** onset near σ ≈ −0.015 (above F ≈ 0.5 pN); melting completes near σ ≈ −1.8.
- **Coexistence is constant-torque** (a Maxwell/lever-rule plateau), the experimental fingerprint of a double well.
- **Hysteresis** is observed after torque/temperature jumps, confirming a finite nucleation barrier (slow kinetics).

**Annotated bibliography.**

- **Strick, Allemand, Bensimon, Bensimon, Croquette 1996, *Science* 271:1835** (10.1126/science.271.5257.1835). Founding magnetic-tweezers assay: a torsionally constrained dsDNA on a magnetic bead; the magnet injects a *counted integer number of turns* n, setting ΔLk = n exactly (the **fixed-Lk ensemble**). Force from transverse Brownian fluctuations via equipartition. "Hat curve" over ±500 turns. *Canonical realization of the ensemble the model must reproduce; winding number n is tracked by construction.*

- **Strick, Allemand, Bensimon, Croquette 1998, *Biophys. J.* 74:2016** (10.1016/S0006-3495(98)77908-1). Above a critical force F_c ≈ 0.5 pN, underwinding stops reducing extension and plateaus because negative turns are absorbed by **local denaturation** rather than writhe — the single-molecule signature of torque melting as two-state coexistence at fixed torque. Overwinding above ~3 pN gives a separate (P-DNA) transition.

- **Allemand, Bensimon, Lavery, Croquette 1998, *PNAS* 95:14152** (10.1073/pnas.95.24.14152). (force, σ) phase diagram. B/denatured coexistence from σ ≈ −0.015 to σ = −1 (a constant-torque line); B/P-DNA coexistence at 3 pN (P-DNA: 2.62 bases/turn, +75% extension). Establishes the melted state as a **distinct phase with its own intrinsic twist** coexisting with B-DNA.

- **Bryant, Stone, Gore, Smith, Cozzarelli, Bustamante 2003, *Nature* 424:338** (10.1038/nature01810). **Rotor-bead tracking (RBT).** A nick defines a free swivel; a rotor bead reports local twist via its orbital angle φ(t), read **continuously / unwrapped** (cumulative rotations, not mod 2π). First *direct* single-molecule torque readout (τ = k_ang·⟨Δφ⟩). Confirms linear twist elasticity, constant-torque melting plateau, and torsional modulus C ≈ 410–460 pN·nm² (twist persistence ~100–110 nm). *Template for twist tracking beyond one turn.*

- **La Porta, Wang 2004, *PRL* 92:190801** (10.1103/PhysRevLett.92.190801). Optical torque wrench: torque from the spin angular momentum of light on a birefringent quartz particle, with continuous angle readout — basis of the constant-torque ensemble.

- **Deufel, Forth, Simmons, Dejgosha, Wang 2007, *Nat. Methods* 4:223** (10.1038/nmeth1013). Angular optical trap (AOT) using a nanofabricated quartz cylinder; torque from optical angular-momentum flux, simultaneous force/extension/angle, fixed-Lk by polarization rotation.

- **Sheinin, Forth, Marko, Wang 2011, *PRL* 107:108102** (10.1103/PhysRevLett.107.108102). **The most quantitatively important paper for the melted basin.** AOT measures F, extension, torque, angle simultaneously while underwinding. Melting torque ≈ −10 pN·nm; melted L-DNA is a *defined left-handed structure* (not an open bubble) with intrinsic twist ≈ −13 bp/turn, bending persistence ≈ 3 nm, torsional modulus ≈ 20 nm, rise 0.48 nm/bp. Strongly sequence-dependent below ~5 pN; sequences converge at high force. *Supplies the explicit soft-melted-basin parameters.*

- **Sheinin, Wang 2009, *PCCP* 11:4800** (10.1039/b901646e). Twist-stretch coupling and torque behavior across the buckling/structural transition; B-DNA torsional stiffness and the extension-twist coupling for the duplex basin.

- **Oberstrass, Fernandes, Bryant 2012, *PNAS* 109:6106** (10.1073/pnas.1113532109). High-resolution RBT on designed sequences resolves **sequence-specific cooperative transitions** (AT-rich strand separation, (GC)ₙ → Z-DNA) as discrete constant-torque plateaus, fit with a statistical-mechanical **Ising/zipper** model: per-bp ΔG_i(seq) + a junction/nucleation penalty (cooperativity = the kinetic barrier) + torsional elastic energy. *Experimental template for sequence-dependent ΔG + barrier.*

- **Oberstrass, Fernandes, Lebel, Bryant 2013, *PRL* 110:178103** (10.1103/PhysRevLett.110.178103). Extracts per-base-pair stability and junction energies; observes "backbending" (negative torsional compliance, a coexistence signature) and breathing dynamics. *Most direct experimental source for per-step ΔG(seq) and nucleation/junction barrier values.*

- **Lebel, Basu, Oberstrass, Tretter, Bryant 2014, *Nat. Methods* 11:456** (10.1038/nmeth.2854). Gold-rotor-bead tracking (AuRBT): >100× faster; maps the free-energy landscape G(angle) = −kBT ln P(angle) from equilibrium fluctuations and **dwell-time/rate (Kramers) kinetics** between basins. Continuous unwrapped twist at high speed.

- **Mosconi, Allemand, Bensimon, Croquette 2009, *PRL* 102:078301** (10.1103/PhysRevLett.102.078301). Magnetic-tweezers torque vs supercoiling; clear buckling/denaturation torque plateau, melting torque ≈ −10 pN·nm; twist via bead angular fluctuations, writhe via Fuller's formula.

- **Gao, Hong, Ye, Inman, Wang 2021, *PRL* 127:028101** (10.1103/PhysRevLett.127.028101). AOT + constant-extension method: twist persistence of extended DNA ≈ 22 nm at very low force, plectonemic ≈ 24 nm. Up-to-date calibration of the duplex torsional stiffness and the linear torque–Lk slope before the plateau. *(Note the apparent tension with the ~100 nm value from Bryant 2003 — these differ in force regime and in whether bend-fluctuation renormalization is included; flagged below.)*

- **Smith, Cui, Bustamante 1996, *Science* 271:795** (10.1126/science.271.5250.795). Original observation of the ~65 pN overstretching plateau (~1.7× elongation) and ssDNA elasticity (persistence ~0.75–1 nm) — the soft melted-state stretching reference (tension, not torque).

---

## 2. Theory & phase behavior under torque / tension (B / denatured / plectoneme coexistence)

The unifying theoretical picture: each structural phase (B, melted L, P, S, plectoneme) is a **separate parabolic free-energy well in twist/linking density plus a creation-energy offset**; the equilibrium is the **convex (Maxwell / common-tangent / lever-rule) lower envelope** at fixed torque or fixed Lk, which produces a **constant-torque coexistence plateau**. This is the macroscopic target the CG-RBP double well must reproduce.

**Annotated bibliography.**

- **Marko, Neukirch 2013, *Phys. Rev. E* 88:062722** (10.1103/PhysRevE.88.062722). **Global force–torque phase diagram.** Per-bp extended free energy of phase i:
  `F_i = [ −g_i + (1/2) kBT C_{f,i} (Δψ − (ψ_i − ψ_B) a_i)² ] a_i + ε_i`
  — a harmonic twisting well centered at the phase's intrinsic helical rotation ψ_i, offset by creation energy ε_i. **Parameters (directly usable):** B: ψ_B = 0.598 rad/bp, a_B = 0.34 nm, C_B = 95 nm, A_B = 45 nm, ε_B = 0. **Melted L:** ψ_L = −2π/16 = −0.393 rad/bp (left-handed), a_L = 0.459 nm, C_L = 19 nm, A_L = 7 nm, ε_L ≈ 2.3–2.5 kBT/bp, melting torque ~10 pN·nm. P-DNA: ψ_P = 2π/3.8, ε_P ≈ 12.5 kBT, torque ~40 pN·nm. Strand-sep energy ε_M = 2.7 kBT + 0.2 kBT·ln(M/150mM). *The single best source of per-phase basin parameters.*

- **Marko, Siggia 1994/1995, *Macromolecules* 27:981; 28:8759; *Phys. Rev. E* 52:2912** (10.1021/ma00130a008; 10.1103/PhysRevE.52.2912). Foundational twistable WLC: bending A ≈ 50 nm, twist C ≈ 75–110 nm, twist-bend coupling; Lk = Tw + Wr; constant-torque two-phase (extended ↔ plectonemic) coexistence; strong-stretching force law. *The single-well elastic baseline.*

- **Marko 2007, *Phys. Rev. E* 76:021926** (10.1103/PhysRevE.76.021926). Two-phase coexistence with Maxwell construction → **constant coexistence torque τ_c**; τ = C_eff·(2π σ/h) in the elastic branch; melting torque τ_m ≈ −10 to −11 pN·nm; per-bp B→melted free energy ε_M ≈ 2.5 kBT/bp. Closest analytic analog to the bimodal torque double-well; relates applied-torque and fixed-Lk ensembles. *Recurs across nearly every slice as the canonical coexistence template.*

- **Neukirch, Marko 2011 (and lineage), extended/plectonemic coexistence.** F = x·Fp(σ_p) + (1−x)·Fs(f, σ_s) with the Lk constraint x σ_p + (1−x) σ_s = σ_total; minimize over partition x → buckling transition, coexistence torque, extension-vs-turns. *Cleanest two-basin fixed-Lk template.*

- **Cocco, Monasson, Marko 1999, *PRL* 83:5178** (arXiv:cond-mat/9904277). **Statistical mechanics of torque-induced denaturation.** Couples Ising-like H-bond opening (m_i ∈ {0,1}) to untwisting; first-order B ↔ denatured d-DNA transition parametrized by σ; torque enters as a Legendre term −Γ·θ tilting the double well; Maxwell construction → bimodal energy and τ_c. Sequence via per-pair ε_i (AT < GC). *The most on-target theory paper for the bimodal torque-melting free energy.*

- **Cocco, Monasson, Marko 2001/2002, *PNAS* 98:8608; *PRE* 65/66** (10.1073/pnas.151257598). Force/torque barriers to unzipping. Per-bp ΔG_AT ≈ 1.1 kBT, ΔG_GC ≈ 3.5 kBT; the duplex↔melted **barrier arises from the dsDNA→ssDNA elastic-stiffness change** (opening pays elastic cost before gaining ss entropy); unzipping force ~12 pN; phase boundary f_u ∝ (Γ − Γ_u)^{1/2}.

- **Bouchiat, Mezard 1998, *PRL* 80:1556** (arXiv:cond-mat/9706050) and 2000, *EPJE* 2:377 (arXiv:cond-mat/9904018). Worm-like-rod-chain in Euler angles; maps the Lk-constrained partition function onto a charged particle in a magnetic-monopole field (writhe = Aharonov-Bohm phase), requiring a short-distance cutoff ≈ helix pitch. Clarifies the Fuller local-writhe formula and its regularization — relevant to stable writhe tracking.

- **Moroz, Nelson 1997/1998, *PNAS* 94:14418; *Macromolecules* 31:6333** (10.1073/pnas.94.26.14418). Entropic elasticity of twist-storing polymers at fixed Lk under tension; **effective torsional persistence C_eff < C** renormalized by bend fluctuations (best fit C ≈ 109 nm); torque saturates at Γ_c at melting; explicit **fixed-Lk ↔ fixed-torque Legendre transform.**

- **Marko & Siggia 1995, *Phys. Rev. E* 52:2912.** Quadratic bend-twist elastic free energy with Lk = Tw + Wr and twist-bend coupling; effective torque τ = C_eff·2πσ/L.

- **Lee, Kornyshev et al. 2011, "Thermomechanics of DNA," arXiv:1101.5182.** Consolidated review of B/L/P competing-phase free energies and parameter ranges.

- **Wereszczynski, Andricioaei 2006, *PNAS* 103:16200** (10.1073/pnas.0603850103). All-atom umbrella sampling + WHAM PMF along an RMSD coordinate between B and P-DNA under combined force + torque; B/P/scP triple point at F ≈ 25.7 pN, τ ≈ 34.8 pN·nm. *Methodological precedent for enhanced-sampling a transition coordinate under an applied-torque ensemble.*

---

## 3. Two-state / Ising / Poland-Scheraga denaturation models and sequence-dependent stability (NN params)

The canonical sequence-dependent superhelical-melting free energy has **three parameters** (Benham): per-bp opening ΔG_i(seq), a junction/nucleation penalty (cooperativity = barrier), and a global quadratic residual-linking term. Poland-Scheraga adds loop entropy; the loop exponent c > 2 (self-avoidance) makes melting first-order/bistable. SantaLucia nearest-neighbor parameters supply ΔG(seq).

**Annotated bibliography.**

- **Benham 1979, *PNAS* 76:3870** (PMID:226985) / **1990, *J. Chem. Phys.* 92:6294** / **1992, *J. Mol. Biol.* 225:835** (10.1016/0022-2836(92)90404-8). Founding two-state superhelical denaturation. Free energy `F = Σ_i b_i x_i + a·(#junctions) + (K/2)(residual linking)²`, with b_i sequence-dependent (~1.0–3.5 kcal/mol per bp), **junction/nucleation a ≈ 10–11 kcal/mol** (the cooperativity/barrier term), soft denatured interstrand torsion, K ≈ 1100·RT/N. *Source of the standard SIDD parameter set and the three-parameter structure that maps onto the fixed-Lk ensemble.*

- **Fye, Benham 1999, *Phys. Rev. E* 59:3408** (10.1103/PhysRevE.59.3408). Exact near-linear-scaling transfer-matrix algorithm that handles the **global quadratic Lk constraint** (separating local two-state DOF from one global twist constraint) — exactly the structure of our fixed-linking-number ensemble. Engine behind WebSIDD.

- **Benham, Bi, Zhang — SIDD/SIST/WebSIDD 2004–2015, *Genome Res.* 14:1575; *Bioinformatics* 31:421** (10.1101/gr.2080004; 10.1101/gr.2575904). Genome-scale per-bp destabilization profiles G(x) at fixed σ (physiological σ ≈ −0.06); demonstrates that sequence-dependent ΔG controls *where* melting nucleates and that finite topological domains are less cooperative than linear DNA.

- **Poland, Scheraga 1966, *J. Chem. Phys.* 45:1456/1464** (10.1063/1.1727786). Loop-entropy two-state model: bound bp weight s = exp(−βΔg_stack); open loop of length l weight Ω(l) ~ μ^l·l^{−c}; cooperativity σ_coop ~ 10^{−5}. c controls transition order. *Statistical-mechanical backbone for treating the open region as a distinct entropic state.*

- **Kafri, Mukamel, Peliti 2000, *PRL* 85:4988** (10.1103/PhysRevLett.85.4988). Self-avoidance gives c ≈ 2.1 > 2 → genuine **first-order** denaturation, justifying two discrete basins with a real barrier rather than a smooth crossover.

- **Vologodskii, Lukashin, Frank-Kamenetskii, Anshelevich 1979, *Biopolymers* 18:1131** (10.1002/bip.1979.360181107); **Vologodskii, Cozzarelli et al. 1992–2017, *J. Mol. Biol.* 227:1224** (10.1016/0022-2836(92)90533-P). Helix-coil + superhelical (fixed-Lk) treatment by Monte Carlo / transfer matrix; Gaussian writhe distribution; supercoiling free energy G(σ) ~ (1100 RT/N)σ²; later versions add explicit bubble/structural transitions and sample nucleation. *Template for sampling the fixed-Lk ensemble and Tw/Wr partition.*

- **Manghi, Palmeri, Destainville 2012, *Phys. Rev. E* / *J. Phys. Condens. Matter* 24:235101** (arXiv:0809.0456). Effective **Ising Hamiltonian** H = −Σ J s_i s_{i+1} − Σ h_i s_i obtained by integrating out chain DOF; sequence-dependent field h_i (from NN ΔG), cooperativity J (junction energy), plus global quadratic superhelical term. *The recipe to reduce a full elastic CG melting model to a phenomenological two-state double well — the conceptual bridge from SE(3) elasticity to an Ising double well.*

- **Wartell, Benight 1985, *Phys. Rep.* 126:67** (10.1016/0370-1573(85)90060-2). Review of Poland-Scheraga helix-coil melting; cooperativity σ_coop (~10^{−5}) and loop exponent c as the physical origin of a tunable nucleation barrier.

- **Cocco, Yan, Léger, Chatenay, Marko 2004, *Phys. Rev. E* 70:011910** (arXiv:cond-mat/0309004). Two-state (paired/separated) model combining NN ΔG with ssDNA elasticity (ss persistence ~1–2 nm) to predict force/torque strand separation. *Shows how to combine per-step NN ΔG with single-stranded elasticity to build the melted basin.*

- **Theodorakopoulos 2019, *Phys. Rev. E* 99:032404** (arXiv:1902.05780). Discrete Kratky-Porod ssDNA; force-induced melting thresholds ~64 pN (free) vs ~111 pN (torsionally constrained), the latter including plectonemic entropy. *Demonstrates the ensemble (free vs fixed-Lk) shifts the transition.*

- **Danilowicz et al. (Prentiss lab) 2003, *PNAS* 100:1694** (arXiv:cond-mat/0310633). Experimental dsDNA↔ssDNA phase boundary F_sep(T) in the temperature–force plane; g(T) = ΔH − TΔS per bp. *Temperature dependence calibration.*

**Sequence-dependent stability parameters (NN):**

- **SantaLucia 1998, *PNAS* 95:1460** (10.1073/pnas.95.4.1460) and **SantaLucia, Hicks 2004, *Annu. Rev. Biophys.* 33:415** (10.1146/annurev.biophys.32.110601.141800). **The canonical NN parameter set.** ΔG(seq) = ΔG_init + Σ ΔG_NN(step) + ΔG_termAT + ΔG_sym. ΔG37 per step (kcal/mol, 1 M NaCl): GC/GC −2.24, CG/CG −2.17, GG/CC −1.84, CA/GT −1.45, GT/CA −1.44, GA/CT −1.30, CT/GA −1.28, AA/TT −1.00, AT/AT −0.88, TA/TA −0.58; initiation +0.98 (G·C end)/+1.03 (A·T end), symmetry +0.43. Stability order GC>CG>GG>GA≈GT≈CA>CT>AA>AT>TA. ΔH/ΔS tables for ΔG(T) = ΔH − TΔS. *(See Breslauer et al. 1986, *PNAS* 83:3746, the historical first NN set, since superseded.)*

- **Owczarzy et al. 2004, *Biochemistry* 43:3537** (10.1021/bi034621r) and **2008, *Biochemistry* 47:5336** (10.1021/bi702363u). Monovalent and Mg²⁺/mixed salt corrections to NN ΔG/Tm. Needed because torque/tweezers buffers are usually not 1 M Na⁺; AT-rich is more salt-sensitive. Mixed-buffer selection by R = √[Mg²⁺]/[Mon].

---

## 4. Coarse-grained DNA models that can melt

**Annotated bibliography.**

*oxDNA family (nucleotide-resolution rigid bodies — same SE(3) philosophy at 1-nt resolution; runs in LAMMPS):*

- **Ouldridge, Louis, Doye 2011, *J. Chem. Phys.* 134:085101** (10.1063/1.3552946). Foundational oxDNA. Each nucleotide is a rigid body (3 sites). U = FENE backbone + LJ excluded volume + orientation-modulated Morse **H-bonding** (its rupture = melting) + 3 stacking terms (nearest, cross, coaxial). **Duplex and melted state are two minima of the SAME potential**, fit so finite-duplex Tm match SantaLucia NN. dsDNA bending ~45–50 nm emerges. *Reference for emergent (not imposed) bistability.*

- **Šulc, Romano, Ouldridge, Rovigatti, Doye, Louis 2012, *J. Chem. Phys.* 137:135101** (10.1063/1.4754132). Sequence dependence: HB strength split AT vs GC, stacking per dinucleotide, fit (VMMC + histogram reweighting) to SantaLucia NN. *The ΔG(seq) recipe for a CG model.*

- **Snodin, Randisi, Mosayebi, Šulc et al. 2015 (oxDNA2), *J. Chem. Phys.* 142:234901** (10.1063/1.4921957). Adds major/minor groove asymmetry (→ nonzero twist-bend coupling, correct supercoiling) and Debye-Hückel salt; twist ~34.1 deg/bp.

- **Matek, Ouldridge, Doye, Louis 2015, *Sci. Rep.* 5:7655** (10.1038/srep07655). **Most directly relevant CG torque-melting study.** Lk = Tw + Wr tracking; σ = (Lk−Lk0)/Lk0 set by clamping ends. Denaturation bubbles have near-zero twist, dumping ~1 turn of undertwist locally; bubbles nucleate in **AT-rich** regions (84–91% AT in tip bubbles) and pin plectonemes. Torque response: linear → overshoot at buckling → plateau (~3 pN·nm under their conditions). Force threshold ~2.5 pN separates plectoneme- from bubble-dominated regimes. Relatively **flat** bubble/plectoneme free-energy landscape (low nucleation barrier, reversible). *(The ~3 pN·nm plateau is much lower than the −10 pN·nm experimental melting torque — a model-specific value; flagged below.)*

- **Matek, Ouldridge, Levy, Doye, Louis 2012, *J. Phys. Chem. B* 116:11616** (10.1021/jp3080755). VMMC + order-parameter (extruded-bp count) sampling of cruciform nucleation under negative supercoiling — methodology for twist-driven structural transitions as nucleation/barrier crossing.

- **Romano, Chakraborty, Doye, Ouldridge, Louis 2013, *J. Chem. Phys.* 138:085101** (10.1063/1.4792252). Overstretching = force-induced melting by unpeeling from free ends (~74 pN model; ~65–68 pN expt), not S-DNA; post-overstretch strand keeps ~57% B-DNA stacking; dF/dT = −0.46 pN/K. *Melted basin under tension.*

- **Henrich, Gutiérrez Fosado, Curk, Ouldridge 2018, *EPJE* 41:57** (10.1140/epje/i2018-11669-8; arXiv:1802.07145). **The LAMMPS CG-DNA implementation of oxDNA.** Documents rigid-body integration via quaternions, conversion to body-frame vectors, the LAMMPS torque convention, and custom integrators `fix nve/dotc/langevin` / `fix nve/dot`; applies torque via virtual traps imposing a target pitch. `seqav` vs `seqdep` switch. *Same codebase — borrow torque handling and the quaternion-Langevin integrator.*

*3SPN family (3 sites/nucleotide; LAMMPS plugin):*

- **Hinckley, Freeman, Whitmer, de Pablo 2013 (3SPN.2), *J. Chem. Phys.* 139:144903** (10.1063/1.4822042). Angle-modulated Morse base pairing; **ε_AT = 16.73, ε_GC = 21.18 kJ/mol**; stacking ~13–15 kJ/mol/step; twist 36 deg/step, rise 3.38 Å. Melting via metadynamics, fit to a logistic two-state curve g(T) = 1/(1+exp[A(T−Tm)]). *Concrete per-pair/per-step ΔG magnitudes and the logistic two-state form.*

- **Freeman, Hinckley, Lequieu, Whitmer, de Pablo 2014 (3SPN.2C), *J. Chem. Phys.* 141:165103** (10.1063/1.4897649). Sequence-dependent ground-state geometry per dinucleotide from X3DNA (Twist, Roll, Tilt, Shift, Slide, Rise) and step-dependent stiffness — the explicit-strand analogue of the RBP intrinsic-S ground state. Runs natively in LAMMPS.

*Peyrard-Bishop(-Dauxois) (PBD) and twist-opening mesoscopic models:*

- **Peyrard, Bishop 1989, *PRL* 62:2755** (10.1103/PhysRevLett.62.2755). Per-bp scalar opening y_n with **Morse on-site V(y) = D(e^{−ay}−1)²** (confining well + flat melted plateau) + harmonic stacking; transfer-integral solution. *Canonical asymmetric bistable on-site form.*

- **Dauxois, Peyrard, Bishop 1993, *Phys. Rev. E* 47:R44/684** (10.1103/PhysRevE.47.R44). Anharmonic stacking W = (k/2)[1 + ρ e^{−α(y_n+y_{n−1})}](y_n−y_{n−1})²: stiffness drops from k(1+ρ) (duplex) to k (melted) once open. **The melted state is much softer and the stiffness change itself drives cooperativity and the barrier.** Typical D ~ 0.03–0.05 eV, a ~ 4.5 Å⁻¹, k ~ 0.025 eV/Å², ρ ~ 1–2, α ~ 0.35 Å⁻¹. *Key lesson: encode the melted basin as a softened-stiffness state, not just an energy offset.*

- **Campa, Giansanti 1998, *Phys. Rev. E* 58:3585** (10.1103/PhysRevE.58.3585). Sequence via D_GC ≈ 1.5·D_AT (3 vs 2 H-bonds).

- **Barbi, Cocco, Peyrard 1999, *Phys. Lett. A* 253:358** (10.1016/S0375-9601(99)00059-6). **Helicoidal twist-opening model.** Two DOF/bp: radial r_n (Morse opening) and **angular φ_n (twist about helix axis)**, coupled by a rigid-backbone geometric constraint. Opening a bubble = a radial **breather locked to a 2π twist kink** in φ. *Because φ is continuous/unbounded, a full-turn untwist is a 2π kink — the conceptual template for storing winding in an angular variable.*

- **Zoli 2020, *Chem. Phys. Lett.* 758:137959** (10.1016/j.cplett.2020.137959). Helicoidal PBD with explicit per-base twist and bending; equilibrium twist set by 10.4 bp/turn; path-integral over (r, θ). Twist-dependent stacking forms; sequence-dependent Morse depths.

- **Choi, Kalosakas, Rasmussen, Bishop et al. 2004, *Nucleic Acids Res.*** (10.1093/nar/gkh1015). Genome-wide PBD with site-dependent D_n, a_n; position-dependent bubble probability/lifetime (AT-rich open first).

- **Weber, Haslam, Essex, Neylon 2009, *Nat. Phys.* / related** (10.1038/nphys1371). Fits PBD Morse depth + stacking to large Tm datasets to produce a **dinucleotide (NN) PBD parameter set** — the bridge from SantaLucia NN ΔG to PBD potentials.

- **Manghi, Destainville 2016, *Phys. Rep.* 631:1** (10.1016/j.physrep.2016.04.001). Authoritative review tying PBD/Morse double-wells, twist-opening coupling, bead-spring rotating-strand mesoscopic models, and Kramers kinetics. *Best single conceptual map for bimodal-energy + barrier-crossing forms.*

*Single-stranded / melted-basin elasticity:*

- **Maffeo, Ngo, Ha, Aksimentiev 2014, *J. Chem. Theory Comput.* 10:2891** (10.1021/ct500193u). ssDNA CG model: persistence length ~0.7–1.5 nm, negligible torsional modulus, FJC/WLC stretching. *Target numbers for the soft melted basin.*

*Twistable WLC / rod elastic continuum (duplex basin):*

- **Skoruppa, Laleman, Nomidis, Carlon 2017, *J. Chem. Phys.* 146:214902** (10.1063/1.4984039). Duplex elastic constants from oxDNA fluctuations mapped to TWLC + Marko-Siggia twist-bend coupling: A ≈ 50 nm, C ≈ 105–118 nm, renormalized twist ~80–110 nm, twist-bend G = 30 nm (oxDNA2) vs ~0 (oxDNA1). Method = covariance/least-squares of (tilt, roll, twist) — directly analogous to the M-matrix harmonic fit.

- **Nomidis, Skoruppa, Carlon, Marko 2019, *Phys. Rev. E* 99:032414** (10.1103/PhysRevE.99.032414). Analytic TWLC + twist-bend coupling: βE = (1/2)[A₁ω₁² + A₂ω₂² + Cω₃² + 2Gω₂ω₃]/a, with Ω = (tilt,roll,twist) = log(R); G renormalizes apparent C. Applied torque ↔ fixed-Lk ensemble.

- **Becker, Everaers 2007, *Phys. Rev. E* 76:021923** (10.1103/PhysRevE.76.021923). Coarse-grains per-step 6×6 stiffness to bend/twist persistence lengths (length-scale-dependent apparent twist stiffness) — the dictionary linking bead-level (10 bp) stiffness to A ~ 50 nm, C ~ 75–110 nm.

- **Brackley, Morozov, Marenduzzo 2014, *J. Chem. Phys.* 140:135103** (arXiv:1404.1704). **The canonical recipe for coding twist into LAMMPS at the bead level.** Each bead is a rigid body with 3 patches → a triad (u_i, f_i, v_i) = exactly T_i. E = (1/2)∫[κ_b(ω₁²+ω₂²) + κ_t ω₃²]ds. LAMMPS "Model 1": bending via an angle style on cos(u_i·u_{i+1}); **twist via TWO dihedral styles** on the patch frames giving the signed discrete twist (α_i+γ_i); a 4th angle aligns u_i with the tangent. Total twist Ω_n = Σ(α_i+γ_i) — **an unwrapped, winding-aware accumulated sum** a single rotation matrix cannot provide. l_p = 50 nm, torsional persistence 60–75 nm. *Implementation blueprint for Hard Problem (1).*

- **Sayar, Avşaroğlu, Kabakçıoğlu 2010, *Phys. Rev. E* 81:041916** (arXiv:0912.0870); **Lionberger/Witz/Stasiak minicircle CG model.** Two-bead-per-step CG with explicit twist; Lk = Tw + Wr enforced/checked, Tw from signed incremental angles (unwrapped), Wr from the Gauss double integral; asymmetric buckling transition. *How to maintain Lk and partition Tw/Wr in a coarse bead model.*

*Topologically constrained mesoscale melting:*

- **Michieletto, Fosado, Marenda (Fosado, Marenda, Michieletto) 2017, *PRL* 119:118002** (10.1103/PhysRevLett.119.118002; arXiv:1703.08367). Large-scale Brownian dynamics at fixed Lk: torque drives **phase coexistence of denatured and intact domains**, with a Landau mean-field theory coupling a denaturation order parameter φ to a supercoiling field σ and a topology-dependent bubble-growth scaling law. The transition is broader/less cooperative than linear melting. *Mesoscale justification for a bistable double well and the fixed-Lk ensemble; untwist is absorbed by melted bubbles.*

---

## 5. Twist / linking-number / writhe representation and CONTINUOUS TRACKING

The field universally tracks twist via the topological decomposition **Lk = Tw + Wr** (White/Cǎlugǎreanu/Fuller), with Lk an integer invariant set by the ensemble and Tw an integrated line integral (never read from a single rotation matrix). This entirely sidesteps the ±180° cap.

**Annotated bibliography.**

- **Cǎlugǎreanu 1961, *Czech. Math. J.* 11:588.** Original differential-geometric statement that self-linking splits into total torsion (twist) + writhe. Cited for provenance.

- **White 1969, *Am. J. Math.* 91:693** (doi:10.2307/2373348). The master identity **Lk = Tw + Wr** for a closed framed ribbon. Tw = (1/2π)∮(u × du/ds)·t ds; Wr = Gauss double integral of the axis with itself. *Any twist lost during melting at fixed Lk must be absorbed by writhe.* (For open clamped segments Lk is not strictly defined → use Fuller's relative writhe + end framing.)

- **Fuller 1971/1978, *PNAS* 68:815; 75:3557** (10.1073/pnas.68.4.815; 10.1073/pnas.75.8.3557). **Reduces the nonlocal writhe double-integral to a single integral relative to a reference curve**, enabling continuous frame-to-frame tracking: Wr = Wr₀ + (1/2π)∫(t₀ × t)·(t₀′ + t′)/(1 + t₀·t) ds. Valid while no antipodal points (t never antiparallel to t₀). For a stretched chain (t₀ = ẑ): dWr = (1/2π)(1 − cos θ)dφ. *The practical engine for tracking writhe of open clamped segments.*

- **Dennis, Hannay 2005, *Proc. R. Soc. A* 461:3245** (arXiv:math-ph/0503012). Reconciles Cǎlugǎreanu/White/Fuller and makes explicit the **integer-jump structure**: twist jumps by 2π when the framing aligns with the principal normal; writhe jumps by ±2 at antipodal points. *Explains exactly why a rotation matrix loses winding number and how to recover a consistent unwrapped twist by smooth framing along arclength and along trajectory time.*

- **Klenin, Langowski 2000, *Biopolymers* 54:307** (10.1002/(SICI)1097-0282). Standard discrete algorithm for **exact writhe of a polygonal chain**: Wr = Σ_{i<j} Ω_ij/2π, with Ω_ij the closed-form signed solid angle of the segment-pair spherical quadrilateral. *The go-to per-frame Wr estimator: feed bead positions, get Wr, check Lk = Tw + Wr each frame.*

- **Levitt 1983, *J. Mol. Biol.* 170:723** (10.1016/S0022-2836(83)80129-6). Earlier segment-pair solid-angle writhe formula — algorithmic ancestor of Klenin-Langowski.

- **Neukirch, Starostin 2008, *Phys. Rev. E* 78:041912** (arXiv:0809.1343). Critically analyzes Fuller's single-integral writhe and its **failure at antipodal points** (plectonemes / tightly writhed configurations); gives corrected formulas and criteria (t·t₀ → −1) for when to fall back to the full Klenin-Langowski double integral. *Tells you when the cheap Fuller estimate breaks.*

- **Bergou, Wardetzky, Robinson, Audoly, Grinspun 2008, *ACM Trans. Graph.* 27:63 (Discrete Elastic Rods)** (10.1145/1360612.1360662). Authoritative algorithm for tracking material-frame twist via **parallel transport (Bishop frame) + a scalar twist angle θ_i**, so total twist = Σθ_i is **unwrapped, continuous, additive** (not capped at π). *The cleanest answer to Hard Problem (1); adopted by many supercoiling codes.*

- **Skoruppa, Carlon 2022, *Phys. Rev. E* 105:044120** (arXiv:2205.07735). Discrete TWLC Monte Carlo at **one bead = 10 bp** (the target resolution). Tw accumulated as the running sum of per-segment twist about the local tangent; Wr from discretized Gauss/Fuller; applied-torque/fixed-Lk by constraining Lk and measuring conjugate torque. *Canonical Tw/Wr/Lk bookkeeping at 10-bp resolution.*

- **Skoruppa, Schiessel et al. 2024, "Multiplectoneme phases," arXiv:2410.23145.** Multiple coexisting plectoneme domains with positional/length entropy and nucleation penalties; lever-rule torque plateau. *Domain-nucleation + entropy bookkeeping transferable to duplex/melted domain coexistence.*

---

## 6. Bistable / double-well free-energy FUNCTIONAL FORMS (the math to borrow)

**Annotated bibliography.**

- **Wiggins, Phillips, Nelson 2005, *Phys. Rev. E* 71:021909** (arXiv:cond-mat/0409003). **Kinkable WLC (KWLC).** Per-joint **two-state sum of Boltzmann factors**: ρ(θ) = exp(−β(ξ/2)θ²) + exp(−βε), i.e. a stiff harmonic basin plus a soft "kinked" basin offset by ε. Kink density ζ = e^{−βε}/(1+e^{−βε}). *The canonical "add two states at the level of Boltzmann weights" route to bistability — directly a stiff duplex + soft melted basin with offset ε ≈ ΔG.*

- **Wiggins et al. 2006, *Nat. Nanotechnol.* 1:137** (10.1038/nnano.2006.63); **Vafabakhsh, Ha 2012, *Science* 337:1097** (10.1126/science.1224139). Experimental fits of the two-state mixture P(θ) = (1−ζ)P_WLC + ζP_kink and cyclization J-factors that *require* a two-state bending model — show how to fit the offset and soft-state stiffness from equilibrium and rate data.

- **Multi-well locally-convex ICNN, arXiv:2506.17242 (2025).** Explicit **log-sum-exp (soft-min)** construction gluing convex basins into one smooth differentiable potential:
  `V(x) = −(1/ρ) log[ (1/N) Σ_i w_i exp(−ρ f_i(x)) ]`, with f_i = (1/2)(Y−Y0_i)ᵀ M_i (Y−Y0_i) + g_i.
  ρ→∞ → min_i f_i (sharp cusp/barrier); finite ρ rounds the crossover and **tunes the barrier height**; w_i / g_i set the asymmetry/ΔG; force = softmax-weighted Σ grad f_i (trivial LAMMPS implementation). *The recommended primary form: two shifted Gaussian wells with a single barrier knob.*

- **Landau-Lifshitz / Ginzburg-Landau quartic.** V(φ) = −(k/2)(φ−φ0)² + (u/4)(φ−φ0)⁴ − h(φ−φ0): minima at ±√(k/u), barrier k²/(4u), asymmetry slope h, curvatures V″ = −k + 3u(φ−φ0)². *Minimal differentiable bistable form with independently tunable barrier, asymmetry, curvature — best for a 1D twist CV and Kramers analysis.* (The two well curvatures are coupled to (k,u); to match independent duplex/melted stiffness use a shifted/piecewise quartic or sextic.)

- **Storm, Nelson 2003, *Phys. Rev. E* 67:051906** (10.1103/PhysRevE.67.051906). B↔S as a **two-state Ising/Zimm-Bragg coexistence**: each segment B (stiff, short) or S (soft, long) with a free-energy difference and a domain-wall (cooperativity) cost = the apparent barrier; transfer-matrix [[e^{−βg_B}, e^{−βJ}],[e^{−βJ}, e^{−βg_S}]]. *Replace extension→twist, B→duplex, S→melted for a direct template.*

- **Argudo, Purohit 2014, *Biophys. J.* 107:2151** (10.1016/j.bpj.2014.09.014). **Quartic double-well + Kramers kinetics** for overstretching. H = [V(u) + Cu]L, V(u) = A₄u⁴ − A₂u² (A₄ ~ 500 pN); minima ±√(A₂/2A₄); linear bias C(F,T) tilts the wells (= applied force/torque conjugate). Kramers escape k = Γ exp(−βE_barrier); reproduces force-step kinetics and asymmetric hysteresis. *The most directly transferable single double-well form; map u→twist, C→applied torque.*

- **Dans, Balaceanu, Pasi, Orozco 2021, *PNAS* 118:e2021263118** (10.1073/pnas.2021263118). Base-pair-step parameters (twist, roll) are frequently **multimodal (BI/BII substates)**; a single Gaussian is insufficient. Multivariate **Ising** model: discrete substates (each a harmonic basin) with sequence-dependent fields h_i (= per-step ΔG) and couplings J_ij. *The bimodal scaffold per step; add the melted basin as a second Ising state.*

- **Dans, Ivani, Orozco 2020, *Nucleic Acids Res.* 48:2099** (10.1093/nar/lkz1227) (and 2023 correlated follow-up). Multi-harmonic CG duplex: per-step energy E = −kBT ln[Σ_k w_k(seq) exp(−(1/2)(X−μ_k)ᵀ K_k (X−μ_k))]. *The engineering realization of a multi-well per-step energy with sequence-dependent weights — i.e. a Gaussian mixture.*

- **Lankaš, Šponer, Langowski, Cheatham 2003, *Biophys. J.* 85:2872** (10.1016/S0006-3495(03)74710-9); **Lankaš et al. 2009 *PCCP* 11:10565** (10.1039/b919565n). Rigid-base-pair harmonic stiffness from MD (K = kBT·covariance⁻¹); documents **bimodal, non-Gaussian** distributions in roll/twist/slide for certain steps (CA/TG, BI/BII) → microscopic justification for a two-Gaussian mixture; group-level ground-state factoring g = s·d (the "Y convention"). *(2009 is also the key methodological reference for marginalizing rigid-base→rigid-base-pair and for why factoring at the group level beats subtracting angles.)*

- **Olson, Gorin, Lu, Hock, Zhurkin 1998, *PNAS* 95:11163** (10.1073/pnas.95.19.11163). Knowledge-based per-dinucleotide intrinsic step parameters (twist ~34–36 deg/step) and 6×6 harmonic force constants — the single-well duplex basin and intrinsic twist Ω₃⁰(seq).

- **Marko 2007, *Phys. Rev. E* 76:021926** (already in §2). Multi-parabola g_phase(σ) = (1/2)C_phase(Tw − Tw0_phase)² + ΔG_phase with Maxwell/convex-hull selection → flat torque plateau; the macroscopic limit of the per-step double well with C_melt ≪ C_ds.

---

## 7. Barrier-crossing / rate / enhanced-sampling methodology and collective variables

**Annotated bibliography.**

*Rate theory:*

- **Kramers 1940, *Physica* 7:284** (10.1016/S0031-8914(40)90098-2). Overdamped escape rate k = (ω_min ω_barrier/2πγ) exp(−βΔF‡), or the mean-first-passage integral k⁻¹ = ∫dx e^{βF(x)}/D(x) ∫dy e^{−βF(y)}. *Converts a computed double-well twist PMF + diffusion coefficient into a melting rate; sets how to tune the barrier for kinetics.*

- **Bell 1978, *Science* 200:618** (10.1126/science.347575). Force(torque)-tilted landscape: k(F) = k₀ exp(βF·x‡); landscape V(x) − F·x. *Applied torque adds −τ·Tw, lowering one barrier; catch-bond extensions give load-switchable asymmetric bistability.*

*Free-energy reconstruction / enhanced sampling:*

- **Torrie, Valleau 1977, *J. Comput. Phys.* 23:187** (10.1016/0021-9991(77)90121-8). Umbrella sampling — bias + reweight, the basis for stratified windows across a melting barrier.

- **Kumar, Rosenberg, Bouzida, Swendsen, Kollman 1992, *J. Comput. Chem.* 13:1011** (10.1002/jcc.540130812). **WHAM.** Optimal unbiasing of umbrella windows → PMF along a CV with error bars to place windows near the barrier.

- **Shirts, Chodera 2008, *J. Chem. Phys.* 129:124105** (10.1063/1.2978177). **MBAR** — binless, asymptotically optimal generalization of WHAM; preferred for multidimensional/continuous CVs (twist + melted-step count) with rigorous uncertainties on the barrier height.

- **Laio, Parrinello 2002, *PNAS* 99:12562** (10.1073/pnas.202427399); well-tempered variant Barducci-Bussi-Parrinello 2008. **Metadynamics** — Gaussian hills fill basins; F(s) ≈ −V_G(s). Ideal for 2D CV (unwrapped twist + melted-step count). Implemented in PLUMED (patches LAMMPS).

- **Dellago, Bolhuis, Csajka, Chandler 1998, *J. Chem. Phys.* 108:1964** (10.1063/1.475562). **Transition path sampling** — harvest reactive trajectories without a good reaction coordinate; committor p_B(x)=0.5 defines the true transition state and **validates whether unwrapped twist alone is a good RC.**

- **Allen, Warren, ten Wolde 2005, *PRL* 94:018104** (10.1103/PhysRevLett.94.018104). **Forward Flux Sampling** — rates and paths from natural dynamics, no good RC needed: k_AB = Φ_{A,0}·Π_i P(λ_{i+1}|λ_i). *The standard for DNA hybridization/melting rates with oxDNA; use a discrete melted-step / native-contact order parameter to sidestep angle periodicity.*

- **Ouldridge, Schreck et al. 2013–2018, oxDNA-FFS.** Concrete FFS on a CG DNA model: order parameter = number of correct base pairs; accelerated diffusion; rates and pathways (zippering, fraying). *Template: replace bp count with melted-step / winding-number CV at 10 bp/bead.*

- **E, Ren, Vanden-Eijnden 2002; Maragliano, Fischer, Vanden-Eijnden, Ciccotti 2006, *J. Chem. Phys.* 125:024106.** **String method in collective variables** — minimum free-energy path, saddle location, committor, and position-dependent free energy + diffusion D(α) for a 1D Smoluchowski rate. *Finds the optimal melting pathway and barrier when twist alone is insufficient.*

- **Chodera, Noé 2014, *Curr. Opin. Struct. Biol.* 25:135** (10.1016/j.sbi.2014.04.002). **Markov state models** + TICA/VAMP — assemble rates from many short sub-barrier trajectories; PCCA+ macrostates (duplex/intermediate/melted); data-driven slow CV discovery.

- **Bonomi, Bussi, Camilloni, Tribello et al. — PLUMED 2009/2019, *CPC* 180:1961; *Nat. Methods* 16:670** (10.1016/j.cpc.2009.05.011). The practical software layer: **LAMMPS `fix plumed`** → RESTRAINT (umbrella), METAD (metadynamics/OPES), reweighting/WHAM, PATHMSD, and **custom CVs via CUSTOM/MATHEVAL** for an unwrapped twist accumulator; bias-exchange replica framework. *Concrete route to apply all methods to the CG-RBP force field without touching the integrator.* **Caveat: PLUMED's built-in TORSION CV is wrapped to (−π,π] — multi-turn twist needs a custom/winding CV.** (`fix colvars` is an alternative engine.)

*Direct barrier-crossing studies on DNA melting (the closest templates):*

- **Sicard, Destainville, Manghi 2015, *J. Chem. Phys.* 142:034903** (arXiv:1405.3867). **The single most on-target paper.** Mesoscopic model: two interacting bead-spring **rotating** strands with a **distance-dependent torsional modulus** κ_φ(ρ) ≈ 450 kBT (duplex) → ~0 (melted) — the two-state stiffness switch. Well-tempered metadynamics with CVs (ρ_max = openness, φ_min = minimal twist; equilibrium twist 0.55 rad) reconstructs the 2D landscape. **Nucleation barrier ΔF_op ≈ 22 kBT, closure ΔF_cl ≈ 13–14 kBT, formation ΔF₀ ≈ 8 kBT for a ~10-bp bubble.** The barrier is of **torsional/elastic origin** — collective untwisting is rate-limiting; κ*_φ(L) = kBT/Var(Φ), Φ = Σ φ_i. **Barrier height scales affinely with β·κ_φ — a directly tunable knob.** Opening times reach ms (Arrhenius/Kramers). *Choose twist as the reaction coordinate; tune the barrier via the bubble torsional modulus.*

- **Dasanna, Destainville, Palmeri, Manghi 2013, *Phys. Rev. E* 87:052703** (10.1103/PhysRevE.87.052703). Brownian dynamics of bubble closure (duplex modulus C ≈ 200–300 kBT, ~0 in bubble); closure limited by a **local Kramers crossing of a torsional barrier** (twist re-winding is the slow coordinate); closure times 0.1–100 µs. *Confirms the winding DOF is the rate-limiting reaction coordinate.*

- **Sicard, Destainville, Rousseau, Tardin, Manghi 2020, *Phys. Rev. E* 101:012403** (10.1103/PhysRevE.101.012403). Supercoiled minicircles at fixed Lk = Tw + Wr; tuning σ and circle size controls opening (ms) and closure (µs–min) times via metadynamics. *Direct fixed-Lk-ensemble template for how the global topological constraint biases the local double well.*

- **van Erp, Cuesta-López, Hagmann, Peyrard 2005, *PRL* 95:218104** (10.1103/PhysRevLett.95.218104). PBD + transition-path/interface sampling: nucleation barrier F(n) vs bubble size n; cooperativity. *Rare-event machinery for bubble barrier crossing.*

- **Liebl, Zacharias 2017, *J. Phys. Chem. B* 121** (10.1021/acs.jpcb.7b07701) (the "Unwinding induced melting" study). All-atom **torque at termini**, RC = unwinding angle / supercoiling density σ (a multi-turn dihedral), umbrella sampling + WHAM. **G(σ) is three-regime: (I) harmonic twist elasticity, (II) abrupt melting transition, (III) flat low-slope post-melting** — the double-basin-with-barrier shape. *The atomistic blueprint for the applied-torque + twist-CV free-energy approach.*

- **"How global DNA unwinding causes non-uniform stress distribution and melting" 2020, PMC7228070.** Same group: G(σ), melting at σ ≈ 0.07, transition cost ~6 kcal/mol; **AT-rich (TATA) melts first, undertwisting to ~21.6 deg/bp vs B-form 34.4**; twist persistence 110.8 nm (AT) / 119.9 nm (GC). *Sequence-dependent localization of the melted bubble.*

- **Neukirch, Marko, Wang lab 2017–2019, *Phys. Rev. E* 95:052401 (torsional buckling kinetics)** (10.1103/PhysRevE.95.052401). Buckling extension-jump as **barrier crossing along a writhe/extension RC**; loop-nucleation barrier; rate insensitive to base-pairing defects while defects pin the buckled state. *Template for casting duplex→melted as Kramers escape over a torque-tilted barrier.*

- **Dittmore, Silver, Neuman 2018, *J. Phys. Chem. B*** (10.1021/acs.jpcb.8b07504). High-speed MT, ΔLk in 0.1-turn steps; two-state **hopping** at buckling; barriers from d(ln τ)/d(Lk); three-minimum landscape (unbuckled/curl/buckled); total barrier ~10 kBT, torque-to-TS ~1.1 pN·nm; mismatches lower buckling Lk but barely change rate (energy/entropy compensation). *How to extract Kramers barriers vs linking number for a torque-driven two-basin transition.*

- **Vlaming/Voorspoels, Skoruppa, Vreede, Carlon 2022/2023 — RBB-NA, *JCTC* 19 / arXiv:2208.10286** (10.1021/acs.jctc.2c00889). PLUMED plugin defining the **12 rigid-base parameters (incl. twist, opening) as CVs** to bias all-atom MD and map the anharmonic G beyond the harmonic well. *Demonstrates rigid-base coordinates are good CVs for biased sampling — directly applicable to biasing an unwrapped-twist / melted-fraction CV.*

- **cgNA+min (Petkevičiūtė, Sharma, Maddocks 2026, *NAR*; arXiv:2411.06036).** Constrained minimization of the cgNA+ Gaussian subject to closure + Lk = Tw + Wr; **torque as a Lagrange multiplier conjugate to integrated twist.** *Template for imposing the fixed-Lk / applied-torque ensemble on a rigid-base elastic energy.*

*Quaternion-Langevin integrators (the LAMMPS engine layer):*

- **Davidchack, Handel, Tretyakov 2009, *J. Chem. Phys.* 130:234101** (10.1063/1.3149788) and **Davidchack, Ouldridge, Tretyakov 2017, *J. Chem. Phys.* 147:224103** (10.1063/1.4999771). Geometric quaternion-Langevin integrator (behind `fix nve/dotc/langevin`): exactly conserves |q|, OU process for angular-momentum noise/friction, anisotropic friction, canonical sampling. *The integrator template for SE(3)/ellipsoid CG-RBP beads. Note: it advances orientation in SO(3) and does NOT store winding — unwrapped twist must be reconstructed externally.*

- **LAMMPS `fix nve/asphere`, `fix addtorque`, `fix addtorque/atom`.** Integrate ellipsoid position/quaternion/angular velocity; apply a (time/variable-dependent) torque to a group's COM or to individual finite-size atoms — **the direct mechanism for a constant applied torque on terminal beads.** Combined with end-clamping (`fix rigid` or harmonic angular springs) → applied-torque or quasi-fixed-Lk ensembles.

---

# A. Approaches to tracking twist beyond ±180° (for the SE(3) CG-RBP / LAMMPS model)

The hard fact (Lavery Curves+ 2009, 10.1093/nar/gkp608; Dennis-Hannay 2005): a rotation matrix gives |Ω| = acos((trR−1)/2) ∈ [0,π], cannot store a winding number (360° twist = identity); the Cayley/rotation-vector (tilt,roll,twist) is valid only for |rotation| < π. Melting untwists by ~1 full turn over a 10-bp bead, so an unwrapped twist is mandatory. Concrete methods, with pros/cons:

**A1. Per-step unwrapped accumulation (recommended primary).** Compute the per-step twist as the third component of the Y-convention deformation Φ_d = log(Sᵀ R) (which already removes the 180° intrinsic-twist branch and keeps the per-step deformation small), then maintain a **running winding counter** along the chain and along the trajectory: add ±2π whenever the bare per-step twist jumps by more than π between successive frames.
- *Pros:* trivial on top of the existing model; at 10-bp resolution the per-step twist change between MD frames is small, so unwrapping is unambiguous; matches the Brackley et al. LAMMPS recipe (Σ(α_i+γ_i) via two dihedral patch styles, naturally unbounded). *Cons:* requires storing per-step state (winding integers) between frames; trajectory-time unwrapping fails if the per-step twist changes by >π in one step (use a small enough dump interval / time step).

**A2. Parallel-transport (Bishop frame) + scalar twist (Discrete Elastic Rods, Bergou 2008).** Carry a zero-twist reference frame by parallel transport bead-to-bead; the material frame is the Bishop frame rotated by accumulated θ_i; total twist = Σθ_i is continuous and additive by construction.
- *Pros:* mathematically clean, the standard in supercoiling codes, no branch-cut handling needed within a frame; integrates naturally with Lk = Tw + Wr. *Cons:* an extra framing layer beyond the triads already carried; the reference-frame transport must be reinitialized consistently each frame.

**A3. Lk = Tw + Wr bookkeeping with discrete writhe (White/Fuller; Klenin-Langowski; Skoruppa-Carlon 2022; minicircle CG).** Compute Wr per frame from bead centerlines (Klenin-Langowski exact solid-angle double sum; Fuller single integral against a straight reference for speed), and recover/cross-check Tw = Lk − Wr at fixed Lk.
- *Pros:* gives an independent topological check each frame; Lk is a natural collective variable; handles writhe absorption of the released turn correctly. *Cons:* Klenin-Langowski is O(N²); Fuller's single integral **fails at antipodal points** (plectonemes/tight writhe; Neukirch-Starostin 2008) and must switch to the full double integral there; Lk is only strictly defined for closed/clamped topology (use Fuller relative writhe + explicit end framing for open clamped segments).

**A4. Discrete-dihedral twist in LAMMPS (Brackley et al. 2014).** Decorate each ellipsoid bead with patches; measure signed twist between consecutive frames with two dihedral styles; accumulate Σ(α_i+γ_i).
- *Pros:* native LAMMPS implementation, signed and unwrapped, already validated against torsional persistence length; directly compatible with `fix nve/asphere`. *Cons:* a dihedral increment can itself wrap if a single step's twist exceeds π (same dump-interval caveat as A1).

**A5. Helicoidal angular order parameter (Barbi-Cocco-Peyrard 1999).** Carry an explicit continuous, unbounded angular variable φ per bead so a melting event is a 2π kink in φ. *Conceptual* — useful as the CV definition; in practice realized by A1/A2.

**Recommendation for the CG-RBP model:** use **A1 (Y-convention per-step deformation + winding counter)** as the working twist, combined with **A3 (Lk = Tw + Wr with Klenin-Langowski Wr)** as a per-frame topology check and as the collective variable for the fixed-Lk ensemble. Treat **A2 (parallel transport)** as the robust fallback if per-step unwrapping proves fragile in tightly writhed configurations. Critically: PLUMED's TORSION CV is wrapped — any biasing must use a custom/winding CV (`CUSTOM`/`MATHEVAL` or a `fix colvars` user variable) fed the unwrapped accumulator, **or** bias the discrete melted-step count instead (Section C), which is periodicity-free.

---

# B. Candidate functional forms for a bimodal duplex/melted free energy

Notation: let `Y` be the per-bead/per-step SE(3) deformation vector (your existing coordinate); `Y0_dup`, `Y0_melt` the basin ground states; `M_dup`, `M_melt` the stiffness matrices (`M_melt ≪ M_dup`); `ΔG(seq)` the sequence-dependent inter-basin offset; `τ` the applied torque; `Tw` the (unwrapped) twist component. **The melted basin's twist component of Y0_melt is shifted by ≈ −2π (one turn) per ~10 bp**, with `M_melt` softened toward ssDNA values (bending persistence ~1–3 nm vs ~50 nm; torsional modulus ~0–20 nm vs ~100 nm).

**B1. Log-sum-exp / soft-min of two Gaussian wells (RECOMMENDED PRIMARY).**
```
E(Y) = −(1/ρ) ln[ exp(−ρ E_dup(Y)) + exp(−ρ (E_melt(Y) + ΔG(seq))) ]
E_dup  = ½ (Y−Y0_dup)ᵀ  M_dup  (Y−Y0_dup)
E_melt = ½ (Y−Y0_melt)ᵀ M_melt (Y−Y0_melt)
```
- **Smooth/differentiable:** analytic gradient = softmax-weighted Σ of basin gradients (cheap LAMMPS force). Source: arXiv:2506.17242; KWLC two-state sum (Wiggins 2005) is the ρ=β special case.
- **Asymmetry:** `ΔG(seq)` raises the melted basin (Section C).
- **Barrier control:** `ρ` tunes the rounding/height at the crossover (ρ→∞ sharp cusp; small ρ merges wells). For **independent** kinetic tuning add an explicit bump `+ B·exp(−(Tw−Tw*)²/2w²)` at the crossover twist.
- **Softened melted elasticity:** enters directly via `M_melt` and the −2π twist shift in `Y0_melt`.

**B2. Asymmetric quartic in a 1D twist CV (for reduced/Kramers analysis).**
```
V(φ) = A₄ φ⁴ + A₂ φ² + A₁ φ,   φ = Tw − Tw0,  A₂ < 0
```
minima ±√(−A₂/2A₄), barrier A₂²/(4A₄), asymmetry A₁ = −τ·(lever arm) (Argudo-Purohit 2014; Landau). *Pros:* fully analytic, ideal for Kramers rate and 1D PMF validation. *Cons:* the two well curvatures are coupled to (A₂,A₄); to match independent duplex vs melted stiffness use a shifted/piecewise quartic or sextic.

**B3. Two-state Ising / Boltzmann-sum per bead (KWLC / Storm-Nelson / Benham).** Give each bead a discrete spin s ∈ {duplex, melted}; partition function sums Boltzmann weights of the two harmonic basins plus a **domain-wall (junction) penalty J** between unlike neighbors = the cooperativity/barrier.
```
H = Σ_i [ s_i E_melt,i + (1−s_i) E_dup,i + s_i ΔG_i(seq) ] + J·(# duplex/melted junctions) − τ·Tw
```
*Pros:* directly encodes cooperativity and the nucleation barrier; equilibrium populations and the torque plateau follow from a transfer matrix; maps onto Benham/SIDD and PLUMED-free melted-step CVs. *Cons:* discrete s needs either a continuous switching coordinate (B4) for MD forces or a hybrid MC/MD scheme.

**B4. Gaussian-mixture / multi-harmonic with a continuous switching coordinate (Dans-Orozco 2020/2021).**
```
E(Y) = −kBT ln[ Σ_k w_k(seq) exp(−½ (Y−μ_k)ᵀ K_k (Y−μ_k)) ]
```
The k=2 case (duplex, melted) is mathematically identical to B1 with w_melt/w_dup = exp(−βΔG); microscopically justified by bimodal MD step distributions (Lankaš; ABC/μABC tetranucleotide BI/BII). *Use this framing to make basin weights sequence-dependent.*

**B5. Morse + anharmonic stacking (PBD; explicit opening coordinate).** V(y) = D(1−e^{−ay})² (harmonic-like duplex well + flat melted plateau), with **state-dependent stacking** W = (k/2)[1 + ρe^{−α(y_n+y_{n−1})}](y_n−y_{n−1})² so the inter-bead stiffness drops k(1+ρ)→k on melting. *Pros:* the canonical microscopic melting form; the stiffness switch generates cooperativity/barrier intrinsically. *Cons:* introduces an extra opening DOF y per bead rather than reusing the twist coordinate; less natural at 10-bp resolution.

**B6. Macroscopic per-phase parabolas + Maxwell construction (Marko 2007; Marko-Neukirch 2013).** g_phase(σ) = (1/2)C_phase(Tw − Tw0_phase)² + ΔG_phase; the convex lower envelope gives the flat coexistence torque. *Not an MD force law* but the **parameter-fixing constraint**: τ_c = (ΔG_melt − ΔG_dup)/(Tw0_dup − Tw0_melt) must equal the experimental −10 to −11 pN·nm. Use to calibrate B1–B4.

**Applied-torque / fixed-Lk in all forms:** add a linear tilt `− τ·Tw` (constant-torque, Gibbs/Legendre) or constrain total Lk = Tw + Wr with a Lagrange multiplier (constant-Lk; cgNA+min). The tilt makes the melted basin the global minimum at the melting torque; coexistence (degenerate basins) defines τ_c.

**Recommendation:** implement **B1** (log-sum-exp of two shifted Gaussians) as the bead/step energy with an explicit barrier bump for independent kinetic tuning; reduce to **B2** in a 1D twist CV for Kramers rates; use **B3** cooperativity (junction penalty) if inter-bead melting cooperativity must be reproduced; fix all parameters against **B6** (τ_c = −10 to −11 pN·nm) and the §1/§2 elastic numbers.

---

# C. Sequence-dependent duplex free energy & melted-state elasticity parameters (concrete numbers/sources)

**C1. Inter-basin offset ΔG(seq) — duplex stabilization per bead.** Sum SantaLucia NN ΔG over the ~9–10 dinucleotide steps in each 10-bp bead:
```
ΔG_bead(seq) = ΔG_init + Σ_{steps in bead} ΔG_NN(step) + ΔG_termAT (+ salt correction)
```
ΔG37 per step (kcal/mol, 1 M Na⁺; SantaLucia 1998 / SantaLucia-Hicks 2004): GC/GC −2.24, CG/CG −2.17, GG/CC −1.84, CA/GT −1.45, GT/CA −1.44, GA/CT −1.30, CT/GA −1.28, AA/TT −1.00, AT/AT −0.88, TA/TA −0.58; initiation +0.98/+1.03, symmetry +0.43. **Per-bead range ≈ 7 kBT (AT-rich) to ~24 kBT (GC-rich)** favoring duplex at zero torque (≈ 0.7–2.4 kcal/mol·bp → ~1.1–3.9 kBT/bp at 300 K). Use ΔG(T) = ΔH − TΔS (full ΔH/ΔS in SantaLucia 1998) at the simulation temperature, then **salt-correct** (Owczarzy 2004 monovalent; Owczarzy 2008 Mg²⁺/mixed; AT-rich is more salt-sensitive). Cross-checks: Cocco-Marko per-bp ΔG_AT ≈ 1.1 kBT, ΔG_GC ≈ 3.5 kBT; Benham per-bp opening 1.0–3.5 kcal/mol; Marko-Neukirch sequence-averaged ε_M ≈ 2.5 kBT/bp with salt term ε_M = 2.7 kBT + 0.2 kBT·ln(M/150mM).

**C2. Nucleation / junction (cooperativity) barrier.** Benham junction a ≈ 10–11 kcal/mol; Poland-Scheraga σ_coop ~ 10⁻⁴–10⁻⁵; **10-bp bubble nucleation barrier ≈ 22 kBT, closure ≈ 13 kBT** (Sicard-Manghi 2015; Manghi-Destainville 2016) — the natural target for the tunable kinetic barrier at 10-bp resolution. Buckling-only barriers are smaller (~10 kBT, Dittmore 2018).

**C3. Duplex basin elasticity (stiff well).** Intrinsic twist ~34.3–36 deg/step (≈10.4–10.5 bp/turn, ψ_B = 0.598 rad/bp); bending persistence A ≈ 45–50 nm; torsional modulus C ≈ 95–110 nm (Marko-Neukirch C_B = 95 nm; Bryant C ≈ 410–460 pN·nm² ≈ 100–110 nm); twist-bend coupling G ≈ 30 nm (oxDNA2/grooved). Sequence-dependent intrinsic shape and stiffness from cgDNA/cgNA+ (Gonzalez-Petkevičiūtė-Maddocks 2013; Sharma-Patelli-Maddocks 2023), Olson 1998, Lankaš 2003, coarse-grained to 10 bp by Skoruppa-Schiessel 2025 (PolyCG, A ≈ 37–53 nm for repeats, correct negative twist-stretch coupling). *The existing model already supplies K(seq), Y0(seq) for this basin.*

**C4. Melted basin elasticity (soft well).** From Sheinin-Wang 2011 / Marko-Neukirch 2013 (L-DNA): **intrinsic twist ≈ −13 bp/turn = −0.393 rad/bp (left-handed; shift ≈ −1 turn per 10 bp)**, bending persistence A_L ≈ 3–7 nm, torsional modulus C_L ≈ 19–20 nm, rise 0.459–0.48 nm/bp. ssDNA-limit alternatives: persistence ~0.7–3 nm, **near-zero torsional modulus** (Maffeo 2014; Sicard-Manghi κ_φ → 0; Cocco et al. ss l_p ~1–2 nm). Sicard-Manghi duplex torsional modulus κ*_φ,ds ≈ 450 kBT (Dasanna 200–300 kBT) drops to ~0 in the bubble — the parameter that sets the barrier height (tunable via β·κ_φ). Practical: set `M_melt` from L-DNA values (twist modulus ~20 nm, bending ~3–7 nm), or push toward ssDNA (~0 twist, ~1–3 nm bending) for a fully denatured bubble; shift the twist component of `Y0_melt` by −2π over the bead.

**C5. Melting torque (calibration target).** τ_melt ≈ −10 to −11 pN·nm (Sheinin-Wang, Bryant, Mosconi, Marko 2007); some single-molecule/elastic-rod analyses report a higher critical magnitude (Vologodskii-Frank-Kamenetskii 2018, τ_c ~ −31 pN·nm). Buckling torque ~+(10–30) pN·nm depending on force/salt. **Flagged discrepancy:** oxDNA (Matek 2015) shows a ~3 pN·nm post-buckling plateau under its specific force/salt — a model-specific value, not the −10 pN·nm thermodynamic melting torque; do not transfer it directly.

---

# Top 15 must-read papers (ranked for this project)

1. **Sicard, Destainville, Manghi 2015, *J. Chem. Phys.* 142:034903** (arXiv:1405.3867) — the closest existing realization of the entire program: bimodal duplex/melted free energy from a distance-dependent torsional modulus, metadynamics with twist as the rate-limiting CV, explicit ~22/13 kBT bubble barriers at 10-bp scale, barrier tunable via β·κ_φ.
2. **Marko, Neukirch 2013, *Phys. Rev. E* 88:062722** — the multi-basin harmonic free energy and every per-phase parameter (B/L/P), melting torque, L-DNA twist −2π/16; the primary source for basin parameters.
3. **Sheinin, Forth, Marko, Wang 2011, *PRL* 107:108102** — direct AOT measurement giving melting torque −10 pN·nm and the full soft melted-basin elasticity (L-DNA twist −13 bp/turn, persistence ~3 nm, torsional ~20 nm).
4. **Marko 2007, *Phys. Rev. E* 76:021926** — two-phase coexistence + Maxwell construction → constant-torque plateau; the analytic template the double well must reproduce and the ensemble link.
5. **Brackley, Morozov, Marenduzzo 2014, *J. Chem. Phys.* 140:135103** (arXiv:1404.1704) — the definitive recipe to encode twist in LAMMPS via patchy rigid-body triads + dihedral styles with an unwrapped twist sum (Hard Problem 1, implementation).
6. **Bergou et al. 2008, *ACM TOG* 27:63 (Discrete Elastic Rods)** — parallel-transport framing for continuous/unwrapped twist beyond ±180° (Hard Problem 1, math).
7. **Cocco, Monasson, Marko 1999, *PRL* 83:5178** (arXiv:cond-mat/9904277) — statistical mechanics of torque-induced denaturation; coupling H-bond opening to untwisting, first-order coexistence, torque-tilted double well.
8. **SantaLucia 1998, *PNAS* 95:1460 / SantaLucia-Hicks 2004** — the NN ΔG(seq) parameters defining the sequence-dependent inter-basin offset.
9. **Skoruppa, Schiessel 2025, *Phys. Rev. Research* 7:013044** (arXiv:2409.05510; PolyCG) — the exact coarse-graining framework underlying the existing model; defines the duplex basin to which the melted basin is added.
10. **Benham 1992, *J. Mol. Biol.* 225:835 + Fye-Benham 1999, *Phys. Rev. E* 59:3408** — the three-parameter superhelical two-state free energy and the exact transfer-matrix algorithm for the global fixed-Lk constraint.
11. **White 1969 + Fuller 1971/1978 + Klenin-Langowski 2000** — Lk = Tw + Wr master identity, single-integral relative writhe, and the discrete per-frame writhe estimator (twist/writhe bookkeeping).
12. **Henrich et al. 2018, *EPJE* 41:57** (arXiv:1802.07145) + **Davidchack-Ouldridge-Tretyakov 2017, *J. Chem. Phys.* 147:224103** — the LAMMPS oxDNA implementation and quaternion-Langevin integrator; the engine layer to reuse.
13. **Allen, Warren, ten Wolde 2005, *PRL* 94:018104 + Laio-Parrinello 2002, *PNAS* 99:12562 (+ Kramers 1940; WHAM/MBAR)** — FFS, metadynamics, and the rate/unbiasing toolkit for getting the double well and melting rates (via PLUMED `fix plumed`).
14. **Matek, Ouldridge, Doye, Louis 2015, *Sci. Rep.* 5:7655** (arXiv:1404.2869) — CG (oxDNA) torque-melting: Lk = Tw + Wr tracking, AT-rich bubble nucleation, plectoneme/bubble competition, fixed-Lk ensemble in practice.
15. **Dans, Balaceanu, Pasi, Orozco 2021, *PNAS* 118:e2021263118 (+ Dans 2020 multimodal)** — multimodal/Ising step energies; the bimodal scaffold with sequence-dependent basin offsets that directly informs the per-step double-well construction.

---

# Gaps / what nobody has done

No published model unifies all the pieces this project requires. (i) The accurate sequence-dependent **rigid-base(-pair) models (cgDNA, cgNA+, MADna, PolyCG, Olson, Lankaš) are strictly single-well, harmonic, and Cayley-vector-bounded** — they have no melted basin, no winding tracking, and explicitly cannot represent twist beyond ±π. (ii) The models that **do melt** are either nucleotide-resolution (oxDNA, 3SPN — far finer than 10 bp/bead, and oxDNA's reported torque plateau ~3 pN·nm does not match the −10 pN·nm thermodynamic melting torque) or two-bead mesoscopic rotating-strand/PBD models that are **not** SE(3) rigid-base-pair force fields and are not sequence-parametrized at the cgDNA level of structural fidelity. (iii) The phenomenological **double-well + barrier with a sequence-dependent ΔG and a softened, −1-turn-shifted melted basin has been written down analytically (Marko, Cocco-Monasson, Manghi) but never implemented as a per-bead 10-bp SE(3) elastic potential in an MD engine** with continuous winding-number tracking. (iv) **Continuous twist tracking beyond one turn inside an SE(3)/quaternion CG-RBP MD model** (as opposed to post-hoc Lk = Tw + Wr analysis of nucleotide models, or graphics-derived parallel transport) appears to be unimplemented for this class of model. (v) **Kramers rates / hysteresis of torque melting at 10-bp resolution with sequence-dependent barriers** exist only at finer (Sicard-Manghi mesoscopic, ~bp) or coarser (Marko thermodynamic) scales — not for a sequence-resolved 10-bp-bead model. The clearest open contribution is therefore exactly the proposed model: a sequence-dependent SE(3) rigid-base-pair double-well at 10 bp/bead, with unwrapped twist (Y-convention deformation + winding counter, cross-checked by Lk = Tw + Wr), a log-sum-exp two-Gaussian energy (NN ΔG offset, L-DNA-soft melted basin shifted by −2π), an independently tunable barrier calibrated to τ_c ≈ −10 pN·nm and ~22 kBT nucleation, run in LAMMPS with quaternion-Langevin integration and applied-torque/fixed-Lk ensembles, with rates from PLUMED metadynamics/FFS + Kramers.

**Caveats / things flagged as uncertain.** Torsional modulus values genuinely differ across sources (C ≈ 95–110 nm from Marko-Neukirch/Bryant vs ~22 nm extended / ~24 nm plectonemic from Gao-Wang 2021) — these reflect different force regimes and whether bend-fluctuation renormalization (C_eff, Moroz-Nelson) is included; pick the value matching your ensemble. The melting-torque magnitude is robustly −10 to −11 pN·nm in tweezers data but a higher critical torque (~−31 pN·nm) appears in some elastic-rod/critical-torque analyses (Vologodskii-Frank-Kamenetskii). The oxDNA ~3 pN·nm plateau is model/condition-specific and should not be used as the melting-torque target. Several identifiers in the source material were incomplete (some Neukirch-Marko, Manghi, and group-attribution entries lacked DOIs, and a few author attributions in the "unwinding-induced melting" all-atom studies were inconsistent across slices — variously attributed to a "Yan group," Liebl-Zacharias, and Mitchell/Harris); these are flagged rather than asserted, and the cited venue/year should be verified before use.
