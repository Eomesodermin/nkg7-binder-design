#!/usr/bin/env python
"""
filter_designs.py — parse Boltz-2 outputs for binder:target complexes,
compute interface metrics, and rank designs.

Metrics (standard de novo binder filters, cf. Bennett et al. 2023, adapted to
Boltz-2 confidence fields which are 0-1 scaled):
  - iptm          : interface pTM (complex/interface confidence; >0.5 pass line)
  - complex_plddt : whole-complex pLDDT 0-1 (fold confidence; >0.7 pass)
  - binder_ptm    : per-chain pTM for the binder chain (fold quality of binder)
  - binder_rmsd   : CA-RMSD between the RFdiffusion backbone and the Boltz-refolded
                    binder (self-consistency; requires the design backbone).

ds>1 AGGREGATION (--aggregate): with diffusion_samples>1 there are N independent
folds per design (model_0..model_{N-1}). Single-fold scores are noisy (measured
ipTM within-5 std ~0.08), so a re-fold/validation pass aggregates the N folds:
  - single       : legacy, model_0 only.
  - mean         : rank on the MEAN of the N per-sample metrics (a lucky single
                   fold no longer wins).
  - mean_penalty : rank on mean_iptm - penalty*iptm_std, penalising designs whose
                   fold is unstable across samples (--stability-penalty, default 1.0).
Aggregated runs emit *_std and n_models columns.

A design PASSES if it clears all configured thresholds. Ranking is by a
composite score (see composite_score) among passers.

Reusable: chain identities and thresholds are CLI args.

Boltz output layout (per input <name>):
  <boltz_out>/boltz_results_<name>/predictions/<name>/
      confidence_<name>_model_0.json   (iptm, complex_plddt, ptm, chains_ptm, ...)
      <name>_model_0.pdb  (or .cif)

Usage:
  python filter_designs.py \
      --boltz-dir 03_boltz --backbone-dir 01_rfdiffusion \
      --binder-chain B --target-chain A \
      --iptm-min 0.5 --plddt-min 0.7 --binder-ptm-min 0.5 --rmsd-max 2.0 \
      --out 04_filtered/ranked.csv
"""
import argparse, json, os, glob, re
import numpy as np


def ca_coords(pdb, chain):
    xs = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
                xs.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xs)


def ca_rmsd(pdb1, chain1, pdb2, chain2):
    """CA-RMSD (Kabsch) between the binder chain in the design backbone and the
    Boltz-refolded complex. Same sequence -> same CA count expected."""
    A, B = ca_coords(pdb1, chain1), ca_coords(pdb2, chain2)
    if len(A) == 0 or len(A) != len(B):
        return np.nan
    A = A - A.mean(0); B = B - B.mean(0)
    U, S, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    A2 = (R @ A.T).T
    return float(np.sqrt(((A2 - B) ** 2).sum(1).mean()))


def binder_ptm(conf, binder_chain_idx):
    cp = conf.get("chains_ptm", {})
    return cp.get(str(binder_chain_idx), None)


def composite_score(row):
    """Higher is better. Boltz fields are 0-1 already."""
    iptm = row["iptm"] or 0
    plddt = row["complex_plddt"] or 0
    bptm = row["binder_ptm"] or 0
    rmsd = row["binder_rmsd"]
    rmsd_term = max(0, 1 - (rmsd / 5.0)) if rmsd == rmsd else 0
    return round(0.45 * iptm + 0.30 * plddt + 0.15 * bptm + 0.10 * rmsd_term, 4)


def find_conf_and_pdb(boltz_dir, name):
    """Locate the Boltz confidence json + model_0 structure for one design."""
    base = os.path.join(boltz_dir, f"boltz_results_{name}", "predictions", name)
    conf = glob.glob(os.path.join(base, "confidence_*_model_0.json"))
    if not conf:  # some layouts drop the boltz_results_ prefix
        conf = glob.glob(os.path.join(boltz_dir, "**", f"confidence_{name}_model_0.json"), recursive=True)
    pdb = glob.glob(os.path.join(base, f"{name}_model_0.pdb")) or           glob.glob(os.path.join(base, f"{name}_model_0.cif"))
    return (conf[0] if conf else None), (pdb[0] if pdb else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz-dir", required=True)
    ap.add_argument("--backbone-dir", default=None,
                    help="RFdiffusion backbone dir for self-consistency RMSD (optional)")
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--binder-chain-idx", type=int, default=1,
                    help="0-based chain index of the binder in the YAML order (target=0, binder=1)")
    ap.add_argument("--iptm-min", type=float, default=0.5)
    ap.add_argument("--plddt-min", type=float, default=0.7)
    ap.add_argument("--binder-ptm-min", type=float, default=0.5)
    ap.add_argument("--rmsd-max", type=float, default=2.0)
    ap.add_argument("--arm", default="")
    ap.add_argument("--aggregate", choices=["single", "mean", "mean_penalty"], default="single",
                    help="ds>1 fold aggregation: single=model_0 only; mean=mean of all samples; "
                         "mean_penalty=mean_iptm - penalty*iptm_std")
    ap.add_argument("--stability-penalty", type=float, default=1.0,
                    help="coefficient on iptm_std subtracted in mean_penalty mode")
    ap.add_argument("--out", required=True)
    ap.add_argument("--extra-boltz-dirs", nargs="*", default=[],
                    help="additional boltz dirs to POOL (full-scale: one per conformer)")
    ap.add_argument("--extra-backbone-dirs", nargs="*", default=[],
                    help="additional RFdiffusion backbone dirs to search for RMSD reference")
    args = ap.parse_args()

    # enumerate designs by their model_0 confidence files (one per design).
    # Full-scale pools across conformers: glob every boltz dir (each recurses into
    # its shard_* subdirs). Design ids carry the arm+conformer tag so they stay unique.
    boltz_dirs = [args.boltz_dir] + list(args.extra_boltz_dirs or [])
    backbone_dirs = ([args.backbone_dir] if args.backbone_dir else []) + list(args.extra_backbone_dirs or [])
    conf_files = []
    for bd in boltz_dirs:
        conf_files += glob.glob(os.path.join(bd, "**", "confidence_*_model_0.json"), recursive=True)
    rows = []
    for cf in sorted(conf_files):
        m = re.search(r"confidence_(.+)_model_0\.json", os.path.basename(cf))
        name = m.group(1) if m else os.path.basename(cf)
        pdir = os.path.dirname(cf)

        # gather ALL diffusion samples for this design (model_0..model_{N-1}).
        # In "single" mode we only use model_0; otherwise we aggregate every model.
        sample_confs = {0: json.load(open(cf))}
        if args.aggregate != "single":
            for mf in glob.glob(os.path.join(pdir, f"confidence_{name}_model_*.json")):
                mm = re.search(r"_model_(\d+)\.json", os.path.basename(mf))
                if mm:
                    idx = int(mm.group(1))
                    if idx not in sample_confs:
                        sample_confs[idx] = json.load(open(mf))
        n_models = len(sample_confs)

        def metric_series(getter):
            vals = [getter(c) for c in sample_confs.values()]
            return np.array([v for v in vals if v is not None], dtype=float)

        iptm_s = metric_series(lambda c: c.get("iptm", c.get("protein_iptm")))
        plddt_s = metric_series(lambda c: c.get("complex_plddt"))
        bptm_s = metric_series(lambda c: binder_ptm(c, args.binder_chain_idx))

        # per-metric point estimate: model_0 in single mode, else the mean.
        def agg(series, model0_val):
            if args.aggregate == "single" or len(series) == 0:
                return model0_val
            return float(np.mean(series))
        conf0 = sample_confs[0]
        iptm = agg(iptm_s, conf0.get("iptm", conf0.get("protein_iptm")))
        plddt = agg(plddt_s, conf0.get("complex_plddt"))
        bptm = agg(bptm_s, binder_ptm(conf0, args.binder_chain_idx))
        iptm_std = float(np.std(iptm_s)) if len(iptm_s) > 1 else 0.0
        plddt_std = float(np.std(plddt_s)) if len(plddt_s) > 1 else 0.0

        # RMSD is measured on the model_0 structure (self-consistency of the backbone
        # is fold-invariant enough; the point estimate stays cheap + comparable).
        pdb = None
        for ext in ("pdb", "cif"):
            hits = glob.glob(os.path.join(pdir, f"{name}_model_0.{ext}"))
            if hits:
                pdb = hits[0]; break
        rmsd = np.nan
        if backbone_dirs and pdb and pdb.endswith(".pdb"):
            # design name is <arm>_<conf>_<backbone#>_s<sample>; backbone file is
            # <arm>_<conf>_<backbone#>.pdb. Search every backbone dir (pooled across
            # conformers) for an exact then suffix-stripped match.
            stem = re.sub(r"_s\d+$", "", name)
            bb = []
            for bd in backbone_dirs:
                bb = (glob.glob(os.path.join(bd, f"{name}.pdb"))
                      or glob.glob(os.path.join(bd, f"{stem}.pdb"))
                      or glob.glob(os.path.join(bd, f"{stem.rsplit('_',1)[0]}.pdb")))
                if bb:
                    break
            if bb:
                rmsd = ca_rmsd(bb[0], args.binder_chain, pdb, args.binder_chain)
        row = dict(
            arm=args.arm, design=name,
            iptm=round(iptm, 4) if iptm is not None else np.nan,
            complex_plddt=round(plddt, 4) if plddt is not None else np.nan,
            binder_ptm=round(bptm, 4) if bptm is not None else np.nan,
            binder_rmsd=round(rmsd, 2) if rmsd == rmsd else np.nan,
            iptm_std=round(iptm_std, 4),
            plddt_std=round(plddt_std, 4),
            n_models=n_models,
            confidence_score=conf0.get("confidence_score"),
            pdb=os.path.basename(pdb) if pdb else "",
        )
        # PASS uses the point estimate (mean in aggregate modes) against thresholds.
        row["pass"] = bool(
            (row["iptm"] >= args.iptm_min if row["iptm"] == row["iptm"] else False)
            and (row["complex_plddt"] >= args.plddt_min if row["complex_plddt"] == row["complex_plddt"] else False)
            and (row["binder_ptm"] >= args.binder_ptm_min if row["binder_ptm"] == row["binder_ptm"] else True)
            and (row["binder_rmsd"] <= args.rmsd_max if row["binder_rmsd"] == row["binder_rmsd"] else True)
        )
        row["composite_score"] = composite_score(row)
        # stability-penalised ranking score: subtract penalty*iptm_std from composite.
        if args.aggregate == "mean_penalty":
            row["rank_score"] = round(row["composite_score"] - args.stability_penalty * iptm_std, 4)
        else:
            row["rank_score"] = row["composite_score"]
        rows.append(row)

    # rank on rank_score (= composite in single/mean modes; composite - penalty*std in mean_penalty)
    rows.sort(key=lambda r: (-int(r["pass"]), -r["rank_score"]))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    import csv
    cols = ["arm", "design", "pass", "rank_score", "composite_score", "iptm", "complex_plddt",
            "binder_ptm", "binder_rmsd", "iptm_std", "plddt_std", "n_models",
            "confidence_score", "pdb"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    npass = sum(r["pass"] for r in rows)
    print(f"[filter] {len(rows)} designs, {npass} pass (aggregate={args.aggregate}, "
          f"mean n_models={np.mean([r['n_models'] for r in rows]):.1f}) -> {args.out}")


if __name__ == "__main__":
    main()
