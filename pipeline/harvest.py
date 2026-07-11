#!/usr/bin/env python
"""
harvest.py — collect validation outputs across all arms of a campaign, merge the
optional consensus / cross-species stages, rank, and build a self-contained
winners bundle.

Runs after submit_campaign.sh finishes (or on whatever has completed so far).
Reusable: reads the same campaign config; nothing target-specific.

Per arm it merges up to three signals, whichever are present on disk:
  03_boltz/          primary Boltz-2 validation           -> iptm/plddt/ptm/rmsd  (always)
  04_chai/           Chai-1 orthogonal consensus          -> chai_iptm, consensus_pass
  05_crossspecies/   fold vs the OTHER species' target    -> cross_iptm, cross_species_pass

Then it copies each arm's winners (passers, else top-N by composite score) into
top_candidates/<arm>/ as complex PDB + binder FASTA + a one-row score CSV, so the
deliverable is self-contained rather than pointers into scratch.

Usage:
  python harvest.py --config configs/nkg7_campaign.yaml --mode smoketest \
      --ws /lustre/.../nkg7_binder_design --out harvest_smoketest [--topn 10]
"""
import argparse, os, glob, subprocess, sys, csv, json, yaml, shutil, re


def binder_seqmap(mpnn_dir):
    """design id (<stem>_s<si>) -> binder sequence, matching stage_boltz naming."""
    m = {}
    for fa in glob.glob(os.path.join(mpnn_dir, "seqs", "*.fa")):
        stem = os.path.basename(fa)[:-3]
        recs = open(fa).read().split(">")
        for si, rec in enumerate(recs[2:]):
            lines = rec.strip().split("\n")
            if len(lines) < 2:
                continue
            seq = "".join(lines[1:]).strip()
            if seq:
                m[f"{stem}_s{si}".replace("/", "_").replace(".", "_")] = seq
    return m


def cross_iptm_map(cross_dir):
    """design id -> ipTM against the cross-species target (from Boltz confidence JSONs)."""
    out = {}
    for cf in glob.glob(os.path.join(cross_dir, "**", "confidence_*_model_0.json"), recursive=True):
        m = re.search(r"confidence_(.+)_model_0\.json", os.path.basename(cf))
        if not m:
            continue
        conf = json.load(open(cf))
        iptm = conf.get("iptm", conf.get("protein_iptm"))
        if iptm is not None:
            out[m.group(1)] = round(float(iptm), 4)
    return out


def chai_map(chai_dir):
    """design id -> {chai_iptm, chai_ptm} from stage_chai's chai_scores.csv."""
    out = {}
    p = os.path.join(chai_dir, "chai_scores.csv")
    if os.path.isfile(p):
        for r in csv.DictReader(open(p)):
            if r.get("chai_iptm") not in ("", None):
                out[r["design"]] = dict(chai_iptm=r["chai_iptm"], chai_ptm=r.get("chai_ptm", ""))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--mode", default="smoketest")
    ap.add_argument("--ws", required=True)
    ap.add_argument("--out", default="harvest")
    ap.add_argument("--topn", type=int, default=10, help="winners to bundle per arm when none pass")
    ap.add_argument("--filter-script", default=os.path.join(os.path.dirname(__file__), "filter_designs.py"))
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    f = cfg["defaults"]["filters"]
    chai_iptm_min = f.get("chai_iptm_min", f["iptm_min"])
    cross_iptm_min = f.get("cross_iptm_min", f["iptm_min"])
    os.makedirs(args.out, exist_ok=True)
    winners_root = os.path.join(args.out, "top_candidates")
    arms = list(cfg["arms"].keys())

    combined = []
    per_arm = []
    for arm in arms:
        a = cfg["arms"][arm]
        rundir = f"{args.ws}/runs/{args.mode}/{arm}"
        boltz_dir = f"{rundir}/03_boltz"
        bb_dir = f"{rundir}/01_rfdiffusion"
        mpnn_dir = f"{rundir}/02_proteinmpnn"
        ranked = os.path.join(args.out, f"{arm}_ranked.csv")
        nconf = len(glob.glob(os.path.join(boltz_dir, "**", "confidence_*_model_0.json"), recursive=True))
        if nconf == 0:
            per_arm.append(dict(arm=arm, n_scored=0, n_pass=0, n_consensus="", n_crossreactive="",
                                best_iptm="", best_score="", status="no boltz output"))
            continue
        subprocess.run([sys.executable, args.filter_script,
               "--boltz-dir", boltz_dir, "--backbone-dir", bb_dir,
               "--binder-chain", a["binder_chain"], "--target-chain", a["target_chain"],
               "--iptm-min", str(f["iptm_min"]), "--plddt-min", str(f["complex_plddt_min"]),
               "--binder-ptm-min", str(f["binder_ptm_min"]), "--rmsd-max", str(f["binder_rmsd_max"]),
               "--arm", arm, "--out", ranked], check=True)
        rows = list(csv.DictReader(open(ranked)))

        # --- merge optional consensus + cross-species signals ---
        chai = chai_map(f"{rundir}/04_chai")
        cross = cross_iptm_map(f"{rundir}/05_crossspecies")
        seqmap = binder_seqmap(mpnn_dir)
        n_consensus = n_cross = 0
        for r in rows:
            d = r["design"]
            r["species"] = a["species"]; r["epitope"] = a["epitope"]
            r["binder_seq"] = seqmap.get(d, "")
            # orthogonal consensus: passes primary AND Chai ipTM over threshold
            ci = chai.get(d, {})
            r["chai_iptm"] = ci.get("chai_iptm", "")
            r["chai_ptm"] = ci.get("chai_ptm", "")
            r["consensus_pass"] = str(bool(
                r["pass"] == "True" and r["chai_iptm"] not in ("", None)
                and float(r["chai_iptm"]) >= chai_iptm_min))
            # cross-species reactivity: native ipTM AND cross ipTM both over threshold
            r["cross_iptm"] = cross.get(d, "")
            native_ok = r["iptm"] not in ("", None) and float(r["iptm"]) >= cross_iptm_min
            r["cross_species_pass"] = str(bool(
                native_ok and r["cross_iptm"] not in ("", None)
                and float(r["cross_iptm"]) >= cross_iptm_min))
            n_consensus += (r["consensus_pass"] == "True")
            n_cross += (r["cross_species_pass"] == "True")
            combined.append(r)

        npass = sum(1 for r in rows if r["pass"] == "True")
        best = rows[0] if rows else {}
        per_arm.append(dict(arm=arm, n_scored=len(rows), n_pass=npass,
                            n_consensus=(n_consensus if chai else ""),
                            n_crossreactive=(n_cross if cross else ""),
                            best_iptm=best.get("iptm", ""), best_score=best.get("composite_score", ""),
                            status="ok"))

        # --- winners bundle for this arm ---
        winners = [r for r in rows if r["pass"] == "True"] or rows[:args.topn]
        if winners:
            wdir = os.path.join(winners_root, arm)
            os.makedirs(wdir, exist_ok=True)
            for r in winners:
                d = r["design"]
                # locate the complex PDB by design id
                hits = glob.glob(os.path.join(boltz_dir, "**", f"{d}_model_0.pdb"), recursive=True)
                if hits:
                    shutil.copy(hits[0], os.path.join(wdir, f"{d}_complex.pdb"))
                if seqmap.get(d):
                    with open(os.path.join(wdir, f"{d}_binder.fasta"), "w") as o:
                        o.write(f">{arm}_{d}_binder\n{seqmap[d]}\n")
            with open(os.path.join(wdir, "scores.csv"), "w", newline="") as fh:
                cols = ["design", "pass", "consensus_pass", "cross_species_pass", "composite_score",
                        "iptm", "complex_plddt", "binder_ptm", "binder_rmsd", "chai_iptm", "cross_iptm", "binder_seq"]
                w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
                w.writeheader(); w.writerows(winners)

    # combined ranked table (passers first, then by composite score)
    combined.sort(key=lambda r: (r.get("pass") != "True", -float(r.get("composite_score") or 0)))
    comb_path = os.path.join(args.out, "all_arms_ranked.csv")
    if combined:
        cols = ["species", "epitope", "arm", "design", "pass", "consensus_pass", "cross_species_pass",
                "composite_score", "iptm", "complex_plddt", "binder_ptm", "binder_rmsd",
                "chai_iptm", "cross_iptm", "confidence_score", "pdb", "binder_seq"]
        with open(comb_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader(); w.writerows(combined)

    sum_path = os.path.join(args.out, "campaign_summary.csv")
    with open(sum_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["arm", "n_scored", "n_pass", "n_consensus",
                                           "n_crossreactive", "best_iptm", "best_score", "status"])
        w.writeheader(); w.writerows(per_arm)

    print(f"[harvest] {len(combined)} designs across {len(arms)} arms")
    print(f"[harvest] per-arm summary -> {sum_path}")
    print(f"[harvest] combined ranked -> {comb_path}")
    print(f"[harvest] winners bundle  -> {winners_root}/<arm>/")
    for r in per_arm:
        extra = ""
        if r["n_consensus"] != "": extra += f" consensus={r['n_consensus']}"
        if r["n_crossreactive"] != "": extra += f" cross-reactive={r['n_crossreactive']}"
        print(f"  {r['arm']:18s} scored={r['n_scored']:3} pass={r['n_pass']:3} "
              f"best_iptm={r['best_iptm']}{extra} ({r['status']})")


if __name__ == "__main__":
    main()
