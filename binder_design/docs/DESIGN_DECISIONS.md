# Design decisions — NKG7 extracellular-loop binder campaign

Every non-obvious methodological choice in this pipeline is recorded here with
its rationale, the alternative that was rejected, and the residual risk it
leaves. This is a living document: add a row when you make a design decision so
the final report can cite it.

---

## 1. Target the extracellular loops only (ECL1, ECL2)

- **Decision.** Design binders exclusively against the two extracellular loops
  of NKG7 — ECL1 (residues 30–60) and ECL2 (residues 113–132) — never the
  transmembrane helices or the cytoplasmic termini.
- **Why.** NKG7 is a 4-TM protein (UniProt Q16617; TM 9–29, 61–81, 92–112,
  133–153). Only the two extracellular loops are accessible to an antibody or
  binder on an intact cell. Every existing NKG7 antibody targets the
  intracellular C-terminus and therefore cannot bind live cells — the whole
  point of this campaign is a **cell-surface-accessible** reagent.
- **Rejected alternative.** Targeting the largest/most-conserved surface
  regardless of topology — rejected because TM and cytoplasmic epitopes are
  biologically inaccessible.
- **Residual risk.** The loops are small (ECL1 31 aa, ECL2 20 aa); interface
  area is the binding challenge (see decision 4).

## 2. Keep the full-length target; restrict the interface to loops via hotspots

- **Decision.** Feed RFdiffusion the **entire** NKG7 model (all 4 TM helices +
  both termini) as the fixed target, and constrain the binder to the loops only
  through `ppi.hotspot_res` (loop residues), rather than excising the loops.
- **Why.** With the TM bundle present, it acts as an implicit **steric wall**:
  RFdiffusion cannot place a binder where the membrane/protein core would be,
  so every design approaches from the physically correct (extracellular) side.
- **Rejected alternative.** Excising ECL1/ECL2 into short peptides and designing
  against those — rejected because a free peptide loses the geometric context
  that holds the loop in place and would permit binders wrapping around faces
  that are buried in the intact protein.
- **Residual risk.** The TM region in the AF model is a static idealization; it
  approximates but does not reproduce the true lipid-embedded environment.

## 3. Design against a loop-conformer ENSEMBLE, not a single AlphaFold snapshot

- **Decision.** Generate a small ensemble of loop conformations (restrained MD:
  TM+cytoplasmic core position-restrained, loops free at 310 K) and distribute
  the design campaign across conformers, preferring binders that score well
  against multiple conformers.
- **Why.** ECL1 and ECL2 carry **no disulfide bonds** — all five cysteines
  (positions 4, 19, 80, 145, 155) are in TM/cytoplasmic regions — so nothing
  staples the loops into a fixed shape; they are conformationally mobile. AF
  gives only one snapshot, and its confidence dips in the loops (per-residue
  pLDDT as low as ~62 in ECL1). A binder cut to complement one frozen loop
  shape may not fit the conformations the loop actually samples.
- **Rejected alternative.** Designing against the single AF model conformation —
  rejected because it over-fits the binder to one uncertain, mobile snapshot.
- **Residual risk.** The ensemble is sampled from a force field starting from a
  *model* (there is no experimental NKG7 structure); it reduces but does not
  eliminate the flexibility risk. Ensemble breadth is a QC output (overlay
  figure) so we can see how much conformational space was actually covered.

## 4. Three epitopes per species: ECL1, ECL2, and composite ECL1+ECL2

- **Decision.** Run three parallel epitope definitions — ECL1 alone, ECL2 alone,
  and a composite patch spanning both loops.
- **Why.** Geometry supports a composite target: the two loops are spatially
  adjacent (centroid–centroid 11.3 Å; closest CA approach 4.6 Å between res 31
  and 127; both project ~+23 Å above the membrane band). The composite patch is
  the largest buildable interface and therefore the highest-feasibility arm; the
  single-loop arms are cheaper hedges and localize the epitope if they succeed.
- **Rejected alternative.** Composite only — rejected because it depends on the
  inter-loop geometry in the AF model being correct; single-loop arms hedge that.
- **Residual risk.** ECL2 alone (20 aa) is the marginal arm — it may not present
  enough area for a robust interface. Kept as a low-cost hedge, flagged as
  lowest-probability.

## 5. Separate human and mouse campaigns; cross-reactivity tested post hoc

- **Decision.** Run the full campaign independently against human NKG7 (Q16617)
  and mouse NKG7 (Q99PA5), then test cross-reactivity in silico (re-fold each
  species' top binders against the other's target).
- **Why.** Loop conservation is moderate — ECL1 61% identical, ECL2 70%
  identical between human and mouse — enough divergence that a human-optimized
  binder will not automatically bind mouse, but with conserved cores
  (`SAHSGLWPT`…`GYIHVTQ` in ECL1; `FFSWS` in ECL2) that make genuinely
  cross-reactive binders plausible. Running both yields species-matched top
  binders **and** flags any that happen to be cross-reactive.
- **Rejected alternative.** Design only against a conserved-consensus epitope to
  force cross-reactivity — rejected because it over-constrains an already small
  interface and would sacrifice the strongest species-specific binders.
- **Residual risk.** Doubles compute. Cross-reactivity is predicted, not
  measured — flagged for wet-lab confirmation.

## 6. No experimental structure — everything is model-based

- **Decision.** Proceed with design against AlphaFold models of NKG7.
- **Why.** No experimental structure of NKG7 exists in the PDB. AF models are
  the best available and carry usable confidence in the loops (mean pLDDT ~85–87).
- **Residual risk.** This is the overarching caveat: target fidelity is bounded
  by AF accuracy for a membrane protein with no template. All downstream
  confidence metrics are conditioned on the model being approximately correct.
  Wet-lab validation is required before any biological conclusion.

## 7. Tooling choices

- **RFdiffusion** for backbone generation with `noise_scale=0` (sharper, more
  designable binder backbones — the recommended binder setting).
- **ProteinMPNN with the soluble-model weights** for sequence design (bias toward
  soluble, well-expressing binders; the target chain is held fixed and only the
  binder chain is redesigned).
- **Boltz-2 (single-sequence co-folding)** for in-silico validation — see
  decision 8 for why Boltz-2 rather than AF2/ColabFold. Single-sequence mode
  (`msa: empty` per chain) is standard for de novo binder filtering and avoids
  MSA-server dependence. Filter on complex pLDDT, interface ipTM, and binder-chain
  pTM, with design-vs-refold CA-RMSD as a self-consistency check.
- **Rejected alternative.** Nanobody/VHH scaffold (RFantibody) as the v1
  modality — deferred; de novo mini-binders have the most mature, highest-yield
  in-silico pipeline and no framework/weights setup. The pipeline's target-prep
  and validation stages are modality-agnostic and can be reused if we add a
  nanobody arm later.

## 8. Validator = Boltz-2 (pip), forced by a marvin hardware constraint

- **Decision.** Use **Boltz-2**, pip-installed into a dedicated conda env, as
  the structure-based validator — not the cluster's ColabFold or Boltz-1
  EasyBuild modules.
- **Why (empirical, not a preference).** marvin's GPU (`mlgpu`, A40) nodes have
  **AMD CPUs with AVX2 but no AVX-512**, while the `easybuild-INTEL` software
  stack is compiled with AVX-512. Both the ColabFold/1.5.2 and Boltz-1/0.4.1
  modules therefore crash instantly with `Illegal instruction` (SIGILL, exit
  132) on these nodes — verified independently for each. Generic pip/conda
  wheels (built for baseline AVX2) run fine, which is why RFdiffusion (pip torch)
  and a pip-installed Boltz-2 both work. Boltz-2 is an AF3-class open-weights
  co-folder that emits the same interface confidence metrics (ipTM, complex
  pLDDT, per-chain pTM) the filter needs.
- **Rejected alternatives.** (a) ColabFold/Boltz-1 modules — non-functional on
  mlgpu (SIGILL). (b) Running AF2 on a different (Intel) GPU partition — not
  pursued because the pip-Boltz-2 route is portable across node types and
  removes the hardware dependency entirely.
- **Residual risk.** Single primary validator. This is now mitigated by an
  orthogonal-consensus stage (Chai-1) — see §9.
- **Reusability note.** On a homogeneous Intel cluster the EasyBuild ColabFold
  module would work and AF2 could be swapped back in; the pipeline isolates the
  validator in one stage script + one parser, so this is a localized change.

## 9. Orthogonal consensus + cross-species + winners bundle as pipeline stages

- **Decision.** Promote three previously-manual post-processing steps into
  built-in, config-toggled pipeline stages that run in the SLURM DAG after the
  primary Boltz-2 validation:
  1. **Orthogonal consensus** (`stage_chai.sh`, `04_chai`) — Chai-1 re-folds the
     same binder:target complexes with an independent all-atom model. A design
     is a *consensus hit* only if it passes Boltz-2 **and** Chai ipTM ≥
     `chai_iptm_min`.
  2. **Cross-species reactivity** (`stage_crossspecies.sh`, `05_crossspecies`) —
     each arm's top-N binders are re-folded against the *other* species' target;
     a design passing ipTM against **both** targets is predicted cross-reactive.
  3. **Winners bundle** (`harvest.py` → `top_candidates/<arm>/`) — the deliverable
     copies each winner's complex PDB + binder FASTA + score row into one folder
     instead of leaving pointers into Lustre scratch.
- **Why.** (1) A single cofolder can score a design highly on model-specific
  artifacts; requiring two independent models to agree is a standard, cheap
  guard for de novo binders. Chai-1 is Apache-2.0, pip-installed (AVX2-safe like
  Boltz — §8), and runs on ESM embeddings without an MSA server, so it needs no
  compute-node egress. (2) The user runs separate human and mouse campaigns and
  explicitly wanted cross-reactivity tested; making it a stage means every
  campaign produces the cross read automatically rather than depending on a
  manual follow-up. (3) A self-contained bundle is what a collaborator or wet-lab
  actually needs.
- **Rejected alternatives.** (a) Leaving these as manual steps — rejected because
  they were easy to forget and the cross-species read was a stated requirement.
  (b) A third validator (AF2/OpenFold3) for consensus — deferred; two-model
  consensus is the usual bar and AF2 is blocked on marvin's AMD nodes anyway (§8).
  (c) Folding *all* designs cross-species — wasteful; only the native top-N are
  worth the cross read, controlled by `crossspecies.topn_*`.
- **Residual risk.** Consensus uses ESM-embedding (no-MSA) Chai, a few ipTM
  points below an MSA-backed run — it is a relative agreement filter, not an
  absolute affinity measure. Cross-species ipTM is a predicted-complex plausibility
  signal, not a binding assay; all cross-reactive calls need experimental
  confirmation. Both stages are `enabled:` flags in the config, so a target where
  they don't apply (single species, no second model available) skips them cleanly.
- **Verification.** Both stages were run end-to-end on marvin A40 nodes for two
  smoke-test arms (human_ecl2, mouse_ecl2): Chai scored all 8 designs/arm,
  cross-species folded 8 complexes/arm against the opposite target, and harvest
  merged `consensus_pass`/`cross_species_pass`/`chai_iptm`/`cross_iptm` into the
  ranked table and wrote the winners bundle. 0 consensus/cross hits at smoke
  scale (T=25), as expected.

## 10. Related work from an adjacent lab (external corroboration + refinement blueprint)

Two 2025–2026 bioRxiv preprints from the Hölzel & Hagelueken groups (Bonn) work
the same problem class (de novo AI minibinders against cancer cell-surface
proteins) with the same toolchain we use, and independently corroborate several
of our design choices. Recorded here for the methods write-up and for the
round-2 refinement discussion.

- **Broske et al. (2025), bioRxiv 2025.05.11.652819** — "Development of AI-designed
  protein binders for detection and targeting of cancer cell surface proteins."
  A discovery + experimental-validation + **benchmarking** workflow (RFdiffusion
  generation with an α-helical-bundle fold filter; mammalian cell-surface + phage
  display screening against PD-L1/CD276/VTCN1). Key result relevant to us:
  **Chai-1 ipTM with ESM embeddings correlates with experimental binding success
  and even predicts the deleterious effect of interface mutations.** This is a
  direct, experimentally-grounded endorsement of our consensus scorer choice
  (§9: Chai-1 + ESM embeddings), from real display data. NOT an evolution method.

- **Tan et al. (2026), bioRxiv 2026.03.04.709551** — "Evolutionary algorithms
  accelerate de novo design of potent Nectin-4-specific cancer biologics."
  A **genetic algorithm (GA) wrapped around the AI design tools** to rescue a
  target (Nectin-4) that plain RFdiffusion under-serves. GA skeleton (select →
  diversify → score, tunable selection stringency) seeded from the best
  RFdiffusion designs, run ~50 generations × 3 independent runs. Two swappable
  mutation operators: **(I) structural — RFdiffusion partial diffusion + ProteinMPNN
  re-sequencing; (II) sequence — direct sequence editing scored by Chai-1 + ESM2.**
  Fitness = adjustable linear combination of binder pTM and complex ipTM; pAE-
  interaction < 5 as the success threshold. **5 GA refinement runs total** (run 1
  seeded from 12 seed designs, run 2 from 104, runs 3–5 refining emerging top
  candidates), ~50 generations each. ~408 A5000 GPU-h produced hundreds of
  designs (92.4% pAE-int < 5); experimentally reached single-digit-nM to sub-nM
  Nectin-4 binders. The GA is explicitly input-agnostic (any upstream designer's
  outputs can seed it).

  **Why this matters for our pipeline.** (a) Their "Option I" mutation operator —
  partial diffusion + ProteinMPNN — IS the round-2 refinement stage we flagged to
  add once winners land; Tan et al. show it works to sub-nM experimentally.
  (b) Their fitness function is Chai-1 + ESM ipTM, which we already run as the
  consensus scorer — so an EvoPro-style loop here should use our existing Chai
  scorer rather than AF2. (c) Their harvested designs seed the GA directly, which
  is exactly how our `harvest/top_candidates/` winners would feed a refinement
  arm. When we revisit refinement, this paper is a validated blueprint for the
  partial-diffusion route (vs. EvoPro's sequence-only, AF2-scored GA).

---

## 11. Pre-flight gates before full-scale: fold-noise + loop-basin convergence

Before committing hundreds of GPU-h to a full-scale (10k+ backbone) campaign,
we ran two cheap pre-flight checks. Both returned consequential results that
reshaped the full-scale design — documenting here so the "why we escalated
before scaling" reasoning is preserved.

### Gate 1 — fold-score sampling noise (is ds=1 ranking trustworthy?)
We re-folded 24 production designs (6 per arm, spanning the composite range) at
Boltz diffusion_samples ds=1 and ds=5, and compared against the original
production ipTM.
- ds=1 vs ds=5 (same fresh run): Spearman rho = 0.72 — sampling noise meaningfully
  shuffles the ranking.
- Production ipTM vs a fresh ds=5 refold: rho = **0.47** — low. Single-fold
  production ipTM is a noisy estimator of a design's "true" score.
- Headline case: `human_ecl2_482_s9`, our best human_ecl2 design, production
  ipTM **0.884** -> ds=5 refold **0.195**. It was largely a lucky fold.
  Conversely `human_ecl1_242_s7` rose from ds=1 0.214 -> ds=5 0.688.
- Within-5 spread: mean per-design ipTM std = 0.076; SEM of the mean falls
  ds1 0.076 -> ds5 0.034 -> ds10 0.024 (curve flattens after ds=5).

**Decision:** the full-scale funnel screens at ds=1 (cheap breadth) then
**re-folds a WIDE shortlist at ds=5**, ranking on the **mean** of the 5 folds
(not best-of-5, which re-introduces the lucky-fold bias), with the per-design
std reported. `filter_designs.py` gained `--aggregate {single,mean,mean_penalty}`
(mean_penalty subtracts penalty*iptm_std from the composite to demote unstable
designs). ds=5 is sufficient for a screening re-fold; higher ds gives
diminishing returns.

### Gate 2 — loop-basin convergence (is one static conformer enough?)
We characterized the flexibility (implicit-MD ensemble) but originally DESIGNED
against a single static AlphaFold conformer. To test whether that is
defensible, we ran 4 independent-seed implicit-MD replicas per species and
compared within-replica loop RMSF against cross-replica basin separation:

| loop | within-rep RMSF | cross-rep max | verdict |
|---|---|---|---|
| human_ECL1 | 0.98 | 2.79 | DIVERGENT |
| human_ECL2 | 0.84 | 1.36 | DIVERGENT |
| mouse_ECL1 | 1.00 | 3.57 | DIVERGENT |
| mouse_ECL2 | 0.82 | 1.09 | borderline |

The seeds do NOT settle into one basin — the loops (ECL1 especially) explore
multiple states. Designing against one static conformer risks
shape-complementarity to a state the loop rarely adopts.

**Decision:** escalate loop characterization before full-scale
(gating rule set in advance: divergence -> escalate). Two independent Tier-3
methods, run in parallel: (3a) explicit-solvent MD (TIP3P/PME/NPT, 50 ns) and
(3b) an orthogonal DL ensemble (AlphaFlow MD+Templates, MSA-backed). Pool all
conformer sources, cluster on loop CA-RMSD to identify basins, and design the
full-scale campaign against **2-3 representative conformers per basin** (chosen
by maximum structural spread), not one snapshot. A basin populated by BOTH
physics and DL is high-confidence. This workflow is captured in the
`protein-loop-flexibility` skill.
