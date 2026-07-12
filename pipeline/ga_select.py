#!/usr/bin/env python
"""
ga_select.py — selection / elitism engine for the partial-diffusion genetic
algorithm (Tan et al. 2026, bioRxiv 2026.03.04.709551, "Option I").

Stdlib only (no numpy/pandas) so it runs in any conda env. Called once per
generation by ga_refine.sbatch. Maintains two persistent CSVs per arm:

  registry.csv   one row per BACKBONE ever produced: backbone_id, pdb_path,
                 best_composite, gen_born, scored. Selection ("who gets mutated
                 next") is over the ENTIRE registry -> global elitism: the k best
                 backbones ever seen are always the parents, so the search can
                 never regress and the seeds are never lost.
  hall.csv       one row per DESIGN (backbone+sequence) ever scored, ranked by
                 composite. This is the final deliverable — the global best
                 designs across all generations.

Two modes:
  --mode init   gen 0: seed the registry from the production winners (their
                known composite is injected; NOT re-folded, saving a whole
                generation of GPU). Selects top-k seeds as gen-1 parents.
  --mode step   gen N>=1: fold results are in --ranked; fold each design back to
                its backbone, update registry + hall, select top-k parents for
                gen N+1, and report whether the global best improved (for the
                early-stop patience check in the driver).

Parent hand-off: selected backbone PDBs are copied into --next-parents-dir with
short stable names g{gen}p{rank}.pdb (names stay bounded across generations),
and their ids are written one-per-line to --seed-list-out for stage_refine.sh.
"""
import argparse, csv, os, re, shutil, sys

def read_csv(path):
    if not path or not os.path.exists(path): return []
    with open(path) as fh: return list(csv.DictReader(fh))

def write_csv(path, rows, cols):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({k: r.get(k, "") for k in cols})

def backbone_of(design):
    """design id <backbone>_s<si>  ->  <backbone>."""
    return re.sub(r"_s\d+$", "", design)

REG_COLS = ["backbone_id", "pdb_path", "best_composite", "best_design", "gen_born", "scored"]
HALL_COLS = ["gen", "arm", "design", "backbone_id", "composite_score", "iptm",
             "complex_plddt", "binder_ptm", "binder_rmsd", "pdb"]

def select_and_stage(registry, k, gen, next_parents_dir, seed_list_out):
    """Pick the k best-composite backbones from the whole registry, copy their
    PDBs into next_parents_dir as g{gen}p{rank}.pdb, write the id list."""
    ranked = sorted(registry, key=lambda r: -float(r["best_composite"]))
    chosen = ranked[:k]
    os.makedirs(next_parents_dir, exist_ok=True)
    ids = []
    for i, r in enumerate(chosen):
        pid = f"g{gen}p{i}"
        dst = os.path.join(next_parents_dir, pid + ".pdb")
        if not os.path.exists(r["pdb_path"]):
            print(f"[ga_select] WARN parent pdb missing: {r['pdb_path']}", file=sys.stderr); continue
        shutil.copyfile(r["pdb_path"], dst)
        ids.append(pid)
    with open(seed_list_out, "w") as fh: fh.write(" ".join(ids) + "\n")
    print(f"[ga_select] staged {len(ids)} parents for gen {gen}: "
          f"best_composite={float(chosen[0]['best_composite']):.4f}")
    return float(chosen[0]["best_composite"])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["init", "step"], required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--gen", type=int, required=True, help="the generation about to be RUN (parents are for this gen)")
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--registry", required=True)
    ap.add_argument("--hall", required=True)
    ap.add_argument("--next-parents-dir", required=True)
    ap.add_argument("--seed-list-out", required=True)
    # init mode
    ap.add_argument("--seed-csv", help="refine_seed_backbones.csv (arm,backbone,...,composite_score)")
    ap.add_argument("--seed-rfd-dir", help="dir holding production seed backbone PDBs <backbone>.pdb")
    # step mode
    ap.add_argument("--ranked", help="this generation's filter_designs CSV")
    ap.add_argument("--backbone-dir", help="dir holding this generation's child backbone PDBs")
    ap.add_argument("--prev-best", type=float, default=-1.0, help="global best composite BEFORE this gen")
    args = ap.parse_args()

    registry = read_csv(args.registry)
    hall = read_csv(args.hall)

    if args.mode == "init":
        seeds = [r for r in read_csv(args.seed_csv) if r["arm"] == args.arm]
        for s in seeds:
            bb = s["backbone"]
            pdb = os.path.join(args.seed_rfd_dir, bb + ".pdb")
            registry.append(dict(backbone_id=bb, pdb_path=pdb,
                                 best_composite=s["composite_score"],
                                 best_design=s.get("design", bb), gen_born="0", scored="1"))
            hall.append(dict(gen="0", arm=args.arm, design=s.get("design", bb), backbone_id=bb,
                             composite_score=s["composite_score"], iptm=s.get("iptm", ""),
                             complex_plddt=s.get("complex_plddt", ""), binder_ptm="",
                             binder_rmsd="", pdb=pdb))
        write_csv(args.registry, registry, REG_COLS)
        write_csv(args.hall, hall, HALL_COLS)
        best = select_and_stage(registry, args.k, args.gen, args.next_parents_dir, args.seed_list_out)
        print(f"GA_BEST={best:.4f}")
        return

    # step mode: ingest this generation's scored designs
    rows = read_csv(args.ranked)
    per_bb = {}   # backbone_id -> best row this gen
    for r in rows:
        comp = r.get("composite_score", "")
        if comp in ("", None): continue
        bb = backbone_of(r["design"])
        c = float(comp)
        if bb not in per_bb or c > per_bb[bb]["_c"]:
            per_bb[bb] = dict(row=r, _c=c)
        hall.append(dict(gen=str(args.gen - 1), arm=args.arm, design=r["design"], backbone_id=bb,
                         composite_score=comp, iptm=r.get("iptm", ""),
                         complex_plddt=r.get("complex_plddt", ""), binder_ptm=r.get("binder_ptm", ""),
                         binder_rmsd=r.get("binder_rmsd", ""),
                         pdb=os.path.join(args.backbone_dir or "", r.get("pdb", ""))))
    # add new backbones to registry
    have = {r["backbone_id"] for r in registry}
    for bb, d in per_bb.items():
        pdb = os.path.join(args.backbone_dir, bb + ".pdb")
        if bb in have:
            for r in registry:
                if r["backbone_id"] == bb and d["_c"] > float(r["best_composite"]):
                    r["best_composite"] = f"{d['_c']:.4f}"; r["best_design"] = d["row"]["design"]
        else:
            registry.append(dict(backbone_id=bb, pdb_path=pdb, best_composite=f"{d['_c']:.4f}",
                                 best_design=d["row"]["design"], gen_born=str(args.gen - 1), scored="1"))
    # keep hall bounded + ranked
    hall = sorted(hall, key=lambda r: -float(r["composite_score"] or 0))[:200]
    write_csv(args.registry, registry, REG_COLS)
    write_csv(args.hall, hall, HALL_COLS)
    best = select_and_stage(registry, args.k, args.gen, args.next_parents_dir, args.seed_list_out)
    improved = best > args.prev_best + 1e-6
    print(f"GA_BEST={best:.4f}")
    print(f"GA_IMPROVED={'1' if improved else '0'}")

if __name__ == "__main__":
    main()
