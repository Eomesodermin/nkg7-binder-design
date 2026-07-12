# Partial-Diffusion Genetic Algorithm (GA) Refinement — `06_refine`

Round-2 refinement stage for the de novo mini-binder pipeline. Reimplements
**"Option I" from Tan et al. 2026** (bioRxiv 2026.03.04.709551, Hölzel &
Hagelueken labs): a genetic algorithm that uses **RFdiffusion partial diffusion
as its mutation operator** to improve winning binders from a first-pass campaign.
Independent reimplementation from the paper's method description — no upstream
code required.

## Why this stage exists

The first-pass NKG7 campaign (6 arms × 500 backbones) produced many
interface-confident designs (767/6000 with Boltz ipTM ≥ 0.5) but almost none
with good **fold quality** (only 12/6000 with complex pLDDT ≥ 0.7). ipTM says
"binds"; pLDDT says "folds well" — the limiter was the fold, not the interface.

A **single** partial-diffusion pass does **not** fix this (calibration, job
26501106: no partial_T beat the seeds). That is expected — one perturbation of an
already-good design usually makes it worse. Tan et al.'s gains come from
**iterating selection over many generations**, which is what this stage does.

## The algorithm

```
gen 0:  seed the registry from the first-pass winners (their known composite is
        injected; NOT re-folded — saves a whole generation of GPU)
repeat for G generations:
  1. SELECT  the top-k backbones from the ENTIRE registry (global elitism)
  2. MUTATE  each: RFdiffusion partial diffusion (diffuser.partial_T) on the
             binder backbone, target held rigid  -> VARIANTS children each
  3. RESEQ   ProteinMPNN re-sequences every child backbone
  4. SCORE   Boltz-2 folds the top-2 seqs/child, filter_designs ranks them
  5. INGEST  update registry (best composite per backbone) + hall-of-fame
  6. STOP    early if the global best hasn't improved for PATIENCE generations
```

### Global elitism — the key correctness property

Selection is over the **entire registry of every backbone ever produced**, seeds
included. The k parents each generation are therefore the k best backbones seen
**so far**, and the final deliverable is the global best across all generations.

**Consequence: the GA output can never be worse than the seeds.** This is the
structural fix for the one-shot regression the calibration exposed, and it is
unit-tested (`ga_select.py`): when a generation's children are all worse than the
current best, the elite parents are retained unchanged and `GA_IMPROVED=0` fires
the early-stop counter.

## Files

| File | Role |
|---|---|
| `stage_refine.sh` | the mutation operator: partial diffusion on a list of seed backbones (`diffuser.partial_T`, target held as motif via contig `[A1-165/0 L-L]`) |
| `ga_select.py` | selection/elitism engine (stdlib only). `--mode init` seeds gen 0 from winners; `--mode step` ingests a generation's scores, updates registry + hall, selects next parents, reports improvement |
| `ga_refine.sbatch` | per-arm SLURM array driver: runs the generation loop (refine → MPNN → Boltz → filter → select) with early stopping |
| `submit_ga.sh` | renders + submits the array (one task/arm, throttled `%2`) + one dependent harvest+email job |
| `ga_harvest.py` | collates per-arm halls into a campaign summary (seed vs final, gain, improved?) + winners bundle |

## Calibrated parameters (NKG7 run)

| Param | Value | Rationale |
|---|---|---|
| `PARTIAL_T` | 10 | calibration: best fold quality without the pt20 pose-drift |
| `K` (parents/gen) | 4 | elitism breadth — enough diversity for selection |
| `VARIANTS` (children/parent) | 8 | 32 child backbones/gen |
| `NUM_SEQ` (MPNN/child) | 8 | sequence diversity per backbone |
| `TOP_SEQS` (Boltz/child) | 2 | fold only the 2 best-scoring seqs — validation, not breadth |
| `DIFF_SAMPLES` | 1 | ds=1 first pass (5× cheaper; refine ds later if needed) |
| `GENERATIONS` | 12 | cap; early-stop usually ends sooner |
| `PATIENCE` | 3 | stop after 3 gens with no global-best improvement |

**Arms:** human_ecl1, human_ecl2, mouse_ecl1, mouse_ecl2. Both composite arms
were dropped after the first pass (weakest on every axis; single-loop epitopes
are the productive targets for NKG7).

## Compute footprint & cluster etiquette

- ~32 child backbones/gen → 64 Boltz folds/gen; ~26 min/generation.
- ~3–5 h wall per arm (early-stop dependent), **~13 GPU-hours total**.
- **2 SLURM submissions total** (1 array + 1 notify) — no per-generation job
  spam. Array throttled `%2` so at most 2 GPUs run at once. Early stopping frees
  GPUs when an arm plateaus. Generations are sequential within a job (gen N+1
  needs gen N's selection), so a single long job per arm is the correct,
  GPU-politest shape rather than many small jobs.
- One completion email via a single `afterany` harvest job (not per-task mail).

## Reusability

Target-agnostic: point `refine_seed_backbones.csv` at any campaign's winners and
set the arm names. `stage_refine.sh` derives the binder length per seed from the
PDB, so no per-target contig edits. `PARTIAL_T` should be re-calibrated per target
(a 2-seed × {5,10,20} sweep, ~45 min) — one partial_T does not transfer blindly.

## Outputs

- `runs/refine_ga/<arm>/registry.csv` — every backbone, best composite, birth gen
- `runs/refine_ga/<arm>/hall.csv` — global-best designs (seeds + all generations)
- `runs/refine_ga/<arm>/<arm>_ga_ranked.csv` — per-arm ranked hall
- `runs/refine_ga/ga_summary.csv` — campaign summary: seed vs final best, gain, improved?
- `runs/refine_ga/ga_top_candidates/<arm>/` — top-N refined complex PDBs
