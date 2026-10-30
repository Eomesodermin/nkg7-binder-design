# Running a campaign

A campaign is defined by one YAML file. The scripts render SLURM job scripts from it and chain
the stages as a dependency graph, so retargeting the pipeline means editing the config, not the
code. Everything below was run on a SLURM cluster with A40 GPUs.

## Cluster-specific settings to adapt

| Setting | Where | Notes |
|---------|-------|-------|
| SLURM account | `SLURM_ACCOUNT` environment variable | Required by the renderers and submit scripts. |
| Partitions | `mlgpu_devel`, `mlgpu_short`, `mlgpu_medium`, `intelsr_short` in `render_jobs.py`, `render_fullscale.py`, `submit_*.sh` and the `.sbatch` templates | These are the names of the cluster this was developed on. |
| Conda activation | `CONDA` constant in `render_jobs.py` and `render_fullscale.py` | A `source .../conda.sh` line for the cluster's Miniforge. |
| Workspace | `--ws` argument, or `$WS` for the standalone `.sbatch` files | Any directory with `targets/ pipeline/ software/ runs/ logs/`. |
| Notifications | `NOTIFY_EMAIL` (optional) | Adds one final harvest-and-email job. |

The standalone `.sbatch` files (`ds_calib`, `md_*`, `refine_calib`) carry a
`YOUR_SLURM_ACCOUNT` placeholder to edit before submitting.

## Environments

Four conda environments back the pipeline (see [`../envs/README.md`](../envs/README.md)):
`SE3nv` (RFdiffusion + ProteinMPNN), `boltz` (Boltz-2), `chai` (Chai-1) and `openmm-env`
(loop ensembles). They are installed with generic AVX2 wheels because the cluster's
AVX-512-compiled module stack crashed (SIGILL) on its AMD GPU nodes
(see [`DESIGN_DECISIONS.md`](DESIGN_DECISIONS.md) §8).

Chai-1 weights (~5 GB) are cached once and pointed to with `CHAI_DOWNLOADS_DIR`, so compute nodes
never re-download them.

## 1. Prepare the target (once per target)

```bash
# on a GPU node, in the openmm-env
python pipeline/make_loop_ensemble.py \
    --pdb targets/NKG7_human_Q16617_AF.pdb \
    --loops 30-60 113-132 \
    --n-conformers 12 --platform CUDA \
    --out targets/ensembles/human
```

This writes `conformer_00..11.pdb`, the loop-conformer ensemble used for the flexibility QC figure
and, optionally, for multi-conformer design. Also write a `<species>_target.seq` file (single-line
amino-acid sequence) in `targets/`, which the Boltz stage uses to build binder:target complexes.

## 2. Edit the config

`configs/nkg7_campaign.yaml` defines the campaign:

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
  rfdiffusion: { num_designs_smoketest: 50, num_designs_production: 500, ... }
  proteinmpnn: { num_seq_per_target_smoketest: 8, num_seq_per_target_production: 16, use_soluble_model: true }
  filters:     { iptm_min: 0.5, complex_plddt_min: 0.7, binder_ptm_min: 0.5, binder_rmsd_max: 2.0 }
```

To target a new protein, add arms with the new `target_pdb`, `hotspot_res` (the epitope residues
the binder should contact) and `contig` (adjust the `A1-N` range to the target length). Nothing
else changes.

## 3. Render the job scripts

```bash
export SLURM_ACCOUNT=my_account
python pipeline/render_jobs.py --config configs/nkg7_campaign.yaml \
    --mode smoketest --ws /path/to/workspace --out $WS/pipeline/rendered
```

`--mode smoketest` uses a small parameter set (50 backbones, 8 sequences, Boltz capped per arm) to
prove the chain end to end in about an hour. `--mode production` uses the full set (500 backbones
per arm, 16 sequences per backbone, top 2 folded). For the multi-conformer full-scale campaign use
`render_fullscale.py` with `pipeline/fullscale_campaign.yaml` and submit it in bounded waves with
`submit_fullscale_waves.sh`.

## 4. Submit

```bash
bash pipeline/submit_campaign.sh $WS smoketest                 # all arms
bash pipeline/submit_campaign.sh $WS production human_ecl1     # a subset
```

Each arm runs as a dependency chain, RFdiffusion → ProteinMPNN → Boltz-2, followed by two
post-validation stages that both depend on Boltz: Chai-1 consensus and cross-species re-fold. Arms
run in parallel, and job IDs are logged to `runs/<mode>/campaign_jobids.tsv`. Either
post-validation stage can be disabled with `consensus.enabled` or `crossspecies.enabled` in the
config.

## 5. Harvest and rank

```bash
python pipeline/harvest.py --config configs/nkg7_campaign.yaml --mode smoketest \
    --ws $WS --out harvest_smoketest --topn 10
```

Outputs:

- `all_arms_ranked.csv`: every design, passers first, with `pass` (Boltz-2), `consensus_pass`
  (Boltz-2 and Chai-1), `cross_species_pass` (native and cross-species ipTM), and the inline
  `binder_seq`.
- `campaign_summary.csv`: per-arm counts and best ipTM.
- `top_candidates/<arm>/`: a self-contained winners bundle (complex PDB, binder FASTA,
  `scores.csv`) that copies the files rather than pointing into scratch.

**Metrics** (Boltz-2 confidence, all 0–1): `iptm` (interface, pass > 0.5), `complex_plddt`
(fold, pass > 0.7), `binder_ptm`, and `binder_rmsd` (design-versus-refold self-consistency, pass
< 2 Å). The composite is `0.45·ipTM + 0.30·pLDDT + 0.15·binder pTM + 0.10·RMSD term`.

## 6. Audit the sequences

Always run the complexity audit before trusting a ranking:

```bash
python pipeline/sequence_complexity.py --designs harvest_production/all_arms_ranked.csv \
    --out-dir harvest_production --plot harvest_production/sequence_complexity.png
```

## Optional: partial-diffusion refinement

`submit_ga.sh` runs a genetic-algorithm refinement of the winners. See
[`GA_REFINEMENT.md`](GA_REFINEMENT.md).
