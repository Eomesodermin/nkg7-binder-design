# Pipeline scripts

| Stage | Scripts |
|-------|---------|
| Target preparation | `make_loop_ensemble.py` (implicit-solvent MD loop conformers), `make_loop_ensemble_explicit.py`, `cluster_basins.py`, `alphaflow/` |
| Job rendering and submission | `render_jobs.py`, `render_fullscale.py`, `submit_campaign.sh`, `submit_fullscale_waves.sh` |
| Design stages (env-var driven) | `stage_rfdiffusion.sh`, `stage_proteinmpnn.sh`, `stage_boltz.sh`, `stage_boltz_sharded.sh`, `stage_chai.sh`, `stage_crossspecies.sh`, `stage_refold.sh` |
| Filtering and ranking | `filter_designs.py`, `harvest.py` |
| Sequence audit | `sequence_complexity.py` |
| Refinement (partial-diffusion GA) | `stage_refine.sh`, `ga_select.py`, `ga_harvest.py`, `ga_refine.sbatch`, `submit_ga.sh` |
| Calibration and MD jobs | `ds_calib.sbatch`, `refine_calib.sbatch`, `md_explicit.sbatch`, `md_replicas.sbatch` |

See [`../docs/running_a_campaign.md`](../docs/running_a_campaign.md) for how the pieces fit
together.
