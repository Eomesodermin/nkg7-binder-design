# De novo mini-binder design pipeline

A config-driven, reusable pipeline for designing de novo mini-binders against a
defined epitope on a target protein, and validating them in silico. Built for
the NKG7 extracellular-loop campaign but **not** NKG7-specific — point it at a
different target PDB + epitope residues and it regenerates the whole campaign.

**Pipeline:** RFdiffusion (backbone generation) → ProteinMPNN (sequence design)
→ Boltz-2 (complex validation) → filter/rank.

Runs on the **marvin** HPC cluster (SLURM, A40 GPUs). See
[`docs/DESIGN_DECISIONS.md`](../docs/DESIGN_DECISIONS.md) for the scientific
rationale behind every methodological choice.

---

## Directory layout

```
binder_design/
  pipeline/
    make_loop_ensemble.py   # (target prep) restrained-MD loop-conformer ensemble
    render_jobs.py          # render per-arm SLURM scripts from a config
    submit_campaign.sh      # submit the chained RFdiff->MPNN->Boltz DAG
    stage_rfdiffusion.sh    # standalone stage runner (env-var driven)
    stage_proteinmpnn.sh    # standalone stage runner
    stage_boltz.sh          # standalone stage runner
    filter_designs.py       # parse Boltz output, compute interface metrics, rank
    harvest.py              # collect all arms -> ranked candidate table
  configs/
    nkg7_campaign.yaml      # the campaign definition (target, epitopes, params)
  docs/
    DESIGN_DECISIONS.md     # scientific rationale (read this first)
  figures/
    loop_ensemble_overlay.png  # loop-flexibility QC figure
```

## The marvin environment (one-time setup, already done for NKG7)

Three conda envs on marvin, each chosen to avoid a hardware pitfall (see
DESIGN_DECISIONS §8):

| env | contents | used by |
|-----|----------|---------|
| `SE3nv` | RFdiffusion 1.1.0 + ProteinMPNN (pip torch 1.12+cu116) | stages 1–2 |
| `boltz` | Boltz-2 2.2.1 (pip torch, CUDA 13) | stage 3 |
| `openmm-env` | OpenMM 8.2 + pdbfixer (conda-forge, CUDA) | loop ensemble |

**Critical hardware note:** marvin's `mlgpu` (A40) nodes have **AMD CPUs with
AVX2 but no AVX-512**, while the `easybuild-INTEL` modules are AVX-512-compiled.
The ColabFold and Boltz-1 **modules SIGILL (crash) on these nodes** — that is why
every ML tool here is **pip-installed** (generic AVX2 wheels) rather than
`module load`ed. RFdiffusion jobs additionally need `module load CUDA/11.7.0`
(provides `libcusparse.so.11` for the cu116 wheels) and `numpy<2`.

Workspace on marvin:
`/lustre/scratch/data/dcorvino_hpc-nkg7_binder_design` (90-day scratch
allocation). Layout: `targets/ pipeline/ software/ runs/ logs/`.

---

## Running a campaign

### 1. Prepare the target (once per target)

```bash
# On a GPU node, in the openmm-env:
python make_loop_ensemble.py \
    --pdb targets/NKG7_human_Q16617_AF.pdb \
    --loops 30-60 113-132 \
    --n-conformers 12 --platform CUDA \
    --out targets/ensembles/human
```

Produces `conformer_00..11.pdb` — the loop-conformer ensemble used for the
flexibility QC figure and (optionally) multi-conformer design. Also write a
`<species>_target.seq` file (single-line amino-acid sequence) in `targets/` —
the Boltz stage reads it to build binder:target complexes.

### 2. Edit the config

`configs/nkg7_campaign.yaml` defines the campaign. Key fields:

```yaml
arms:
  human_composite:
    species: human
    epitope: composite
    target_pdb: targets/NKG7_human_Q16617_AF.pdb
    target_chain: A          # chain held fixed
    binder_chain: B          # RFdiffusion writes the binder as chain B
    hotspot_res: [A36, A45, A47, A50, A117, A122, A128, A131]  # epitope residues
    contig: "[A1-165/0 50-65]"   # target fixed (1-165) + binder length range
defaults:
  rfdiffusion: { num_designs_smoketest: 50, num_designs_production: 2000, ... }
  proteinmpnn: { num_seq_per_target_smoketest: 8, ..., use_soluble_model: true }
  validator:   { tool: boltz2, diffusion_samples_smoketest: 1, ... }
  filters:     { iptm_min: 0.5, complex_plddt_min: 0.7, binder_ptm_min: 0.5, binder_rmsd_max: 2.0 }
```

**To target a new protein:** add arms with the new `target_pdb`, `hotspot_res`
(the epitope residues you want the binder to contact), and `contig` (adjust the
`A1-N` range to the target length). Nothing else changes.

### 3. Render the SLURM scripts

```bash
# in the SE3nv env (has pyyaml)
python render_jobs.py \
    --config configs/nkg7_campaign.yaml \
    --mode smoketest \                # or: production
    --ws /lustre/scratch/data/dcorvino_hpc-nkg7_binder_design \
    --out $WS/pipeline/rendered
```

`--mode smoketest` uses the small parameter set (50 backbones, 8 seqs, Boltz
capped to 8 complexes/arm) to prove the chain end-to-end in ~1 h.
`--mode production` uses the full set (2000 backbones, 16 seqs, no Boltz cap).

### 4. Submit the campaign

```bash
bash submit_campaign.sh $WS smoketest          # all arms
bash submit_campaign.sh $WS production human_composite mouse_composite  # subset
```

Each arm runs as a dependency chain: RFdiffusion → (afterok) ProteinMPNN →
(afterok) Boltz-2. Arms run in parallel. Job IDs are logged to
`runs/<mode>/campaign_jobids.tsv`. Watch with `squeue -u $USER`.

### 5. Harvest & rank

```bash
python harvest.py --config configs/nkg7_campaign.yaml --mode smoketest \
    --ws $WS --out harvest_smoketest
```

Produces:
- `all_arms_ranked.csv` — every design across all arms, passers first, ranked by
  composite score.
- `campaign_summary.csv` — per-arm counts (scored, passing, best ipTM).
- `<arm>_ranked.csv` — per-arm detail.

**Metrics** (Boltz-2 confidence, all 0–1): `iptm` (interface, pass > 0.5),
`complex_plddt` (fold, pass > 0.7), `binder_ptm`, `binder_rmsd` (design-vs-refold
self-consistency, pass < 2 Å). `composite_score = 0.45·iptm + 0.30·plddt +
0.15·binder_ptm + 0.10·rmsd_term`.

---

## Scaling from smoke test to production

The smoke test (`--mode smoketest`) proves every stage runs and produces
correctly-shaped outputs; the design metrics are not expected to be good (50
backbones at T=25 gives low-complexity binders). For a real campaign:

- Use `--mode production` (2000 backbones/arm, T=50, 16 seqs, all complexes
  validated).
- Move RFdiffusion/Boltz jobs to `mlgpu_short` (8 h) or `mlgpu_medium` (24 h)
  — the render sets this automatically for `production`.
- Expect the interface filter to pass a small fraction; the composite arm is
  the highest-probability, ECL2-alone the lowest (DESIGN_DECISIONS §4).
- Follow up passing designs with a second orthogonal validator (e.g. Chai-1)
  for consensus, and cross-species re-folding for cross-reactivity.

## Cross-species reactivity check

After harvesting, re-fold each species' top binders against the *other*
species' target (swap the `target` sequence in the Boltz YAML). A binder that
keeps ipTM > 0.5 against both is predicted cross-reactive. Loop conservation is
moderate (ECL1 61 %, ECL2 70 % human↔mouse), so some cross-reactive hits are
plausible but must be confirmed experimentally.
