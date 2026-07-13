#!/usr/bin/env python
"""
cluster_basins.py — pool loop-conformer ensembles from multiple sources,
identify conformational BASINS by clustering on loop CA-RMSD, and max-spread
downselect a small set of design representatives spanning the basins.

Runs after the characterisation jobs finish. Consumes conformer PDBs from any
mix of sources (implicit MD, explicit-solvent MD, AlphaFlow DL) — each source is
a directory of PDBs sharing the target's residue numbering. All conformers are
superposed on the fixed CORE (residues OUTSIDE the loop windows), then compared
by loop-only CA-RMSD, so clustering reflects loop movement, not rigid drift.

BASIN IDENTIFICATION: hierarchical clustering (average linkage) on the loop-RMSD
matrix, cut at --basin-cutoff Angstrom. Each cluster = one basin. The per-basin
source composition is reported: a basin populated by BOTH MD and DL is high-
confidence (independent methods agree it exists).

DOWNSELECTION (--n-select): pick representatives by MAX STRUCTURAL SPREAD, not
arbitrary snapshots. Greedy farthest-point sampling on the loop-RMSD matrix seeded
at the medoid of the largest basin, then repeatedly add the conformer maximally
distant from those already chosen. This guarantees the chosen designs span the
observed conformational range. Optionally forces >=1 representative per basin
(--per-basin) so no distinct state is missed.

Reusable: loop windows, sources, cutoff, and selection count are all CLI args.

Usage:
  python cluster_basins.py --loops 30-60 113-132 \
      --source implicit:runs/md_replicas/human_rep1 ... \
      --source explicit:runs/md_explicit/human \
      --source dl:runs/alphaflow/human \
      --n-select 8 --per-basin --out runs/basins/human
"""
import argparse, glob, os, json
import numpy as np


def parse_windows(specs):
    return [(int(s.split("-")[0]), int(s.split("-")[1])) for s in specs]


def in_windows(resid, windows):
    return any(a <= resid <= b for a, b in windows)


def read_ca(pdb):
    """Return {resid: xyz} for CA atoms (first model / chain A-agnostic)."""
    ca = {}
    for ln in open(pdb):
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            try:
                resid = int(ln[22:26])
            except ValueError:
                continue
            ca[resid] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        elif ln.startswith("ENDMDL"):
            break
    return ca


def kabsch(P, Q):
    """Rotation aligning P onto Q (both centered)."""
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loops", nargs="+", required=True)
    ap.add_argument("--source", action="append", required=True,
                    help="LABEL:DIR  (e.g. explicit:runs/md_explicit/human); repeatable")
    ap.add_argument("--basin-cutoff", type=float, default=2.0,
                    help="loop CA-RMSD (A) linkage cutoff defining separate basins")
    ap.add_argument("--n-select", type=int, default=8)
    ap.add_argument("--per-basin", action="store_true",
                    help="force at least one representative from every basin")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    windows = parse_windows(args.loops)
    os.makedirs(args.out, exist_ok=True)

    # --- load all conformers from all sources ---
    confs = []   # list of dict(path, label, ca)
    for spec in args.source:
        label, d = spec.split(":", 1)
        for pdb in sorted(glob.glob(os.path.join(d, "*.pdb"))):
            if os.path.basename(pdb) == "prepared.pdb":
                continue
            ca = read_ca(pdb)
            if ca:
                confs.append(dict(path=pdb, label=label, ca=ca))
    if len(confs) < 2:
        raise SystemExit(f"need >=2 conformers, found {len(confs)}")
    print(f"[basins] loaded {len(confs)} conformers from {len(args.source)} sources", flush=True)

    # residues common to all conformers
    common = set(confs[0]["ca"])
    for c in confs[1:]:
        common &= set(c["ca"])
    core = sorted(r for r in common if not in_windows(r, windows))
    loop = sorted(r for r in common if in_windows(r, windows))
    print(f"[basins] {len(core)} core residues (superpose), {len(loop)} loop residues (compare)", flush=True)

    # --- superpose every conformer on the core of conformer 0, collect loop coords ---
    ref_core = np.array([confs[0]["ca"][r] for r in core]); ref_core_c = ref_core - ref_core.mean(0)
    ref_core_mean = ref_core.mean(0)
    loopmat = []
    for c in confs:
        cc = np.array([c["ca"][r] for r in core]); ccm = cc.mean(0)
        R = kabsch(cc - ccm, ref_core_c)
        L = np.array([c["ca"][r] for r in loop])
        Lt = (L - ccm) @ R.T + ref_core_mean
        loopmat.append(Lt.ravel())
    loopmat = np.array(loopmat)

    # --- pairwise loop CA-RMSD matrix ---
    n = len(confs); nres = len(loop)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.sqrt(((loopmat[i] - loopmat[j]) ** 2).reshape(nres, 3).sum(1).mean())
            D[i, j] = D[j, i] = d

    # --- hierarchical clustering into basins ---
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    Z = linkage(squareform(D, checks=False), method="average")
    basins = fcluster(Z, t=args.basin_cutoff, criterion="distance")
    nb = len(set(basins))
    print(f"[basins] {nb} basins at cutoff {args.basin_cutoff} A", flush=True)

    # per-basin composition (which methods populate it)
    basin_info = {}
    for b in sorted(set(basins)):
        idx = [i for i in range(n) if basins[i] == b]
        labels = {}
        for i in idx:
            labels[confs[i]["label"]] = labels.get(confs[i]["label"], 0) + 1
        # medoid = conformer with min mean intra-basin distance
        medoid = min(idx, key=lambda i: D[i][idx].mean())
        basin_info[int(b)] = dict(size=len(idx), sources=labels,
                                  medoid=confs[medoid]["path"], members=idx)
        print(f"  basin {b}: n={len(idx)} sources={labels} medoid={os.path.basename(confs[medoid]['path'])}")

    # --- max-spread downselection (greedy farthest-point) ---
    # seed at medoid of largest basin
    largest = max(basin_info.values(), key=lambda v: v["size"])
    seed = min(largest["members"], key=lambda i: D[i][largest["members"]].mean())
    chosen = [seed]
    if args.per_basin:  # ensure each basin has >=1 rep (its medoid)
        for b, info in basin_info.items():
            med = min(info["members"], key=lambda i: D[i][info["members"]].mean())
            if med not in chosen:
                chosen.append(med)
    while len(chosen) < args.n_select and len(chosen) < n:
        # add the conformer maximally distant from the current chosen set
        rest = [i for i in range(n) if i not in chosen]
        nxt = max(rest, key=lambda i: min(D[i][c] for c in chosen))
        chosen.append(nxt)
    chosen = chosen[:max(args.n_select, len(set(basins)) if args.per_basin else args.n_select)]

    # --- write outputs ---
    sel = []
    for rank, i in enumerate(chosen):
        sel.append(dict(rank=rank, conformer=os.path.basename(confs[i]["path"]),
                        source=confs[i]["label"], basin=int(basins[i]), path=confs[i]["path"]))
    summary = dict(n_conformers=n, n_basins=nb, basin_cutoff=args.basin_cutoff,
                   loop_residues=loop, core_residues=len(core),
                   basins={str(k): dict(size=v["size"], sources=v["sources"],
                                        medoid=os.path.basename(v["medoid"]))
                           for k, v in basin_info.items()},
                   selected=sel)
    json.dump(summary, open(os.path.join(args.out, "basin_summary.json"), "w"), indent=2)
    np.save(os.path.join(args.out, "loop_rmsd_matrix.npy"), D)
    with open(os.path.join(args.out, "conformer_index.tsv"), "w") as fh:
        fh.write("idx\tlabel\tbasin\tpath\n")
        for i, c in enumerate(confs):
            fh.write(f"{i}\t{c['label']}\t{basins[i]}\t{c['path']}\n")
    # copy chosen representatives into out/representatives/
    repdir = os.path.join(args.out, "representatives"); os.makedirs(repdir, exist_ok=True)
    import shutil
    for s in sel:
        shutil.copyfile(s["path"], os.path.join(repdir, f"rep{s['rank']:02d}_basin{s['basin']}_{s['source']}.pdb"))
    print(f"[basins] selected {len(sel)} representatives across {nb} basins -> {repdir}", flush=True)
    for s in sel:
        print(f"  rep{s['rank']:02d}  basin {s['basin']}  {s['source']}  {s['conformer']}")
    print("CLUSTER_BASINS_DONE", flush=True)


if __name__ == "__main__":
    main()
