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
