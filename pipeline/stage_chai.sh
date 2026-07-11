#!/bin/bash
# stage_chai.sh — Chai-1 orthogonal-consensus refold of designed binder:target
# complexes (one arm).
#
# Runs the SAME binder:target complexes as the Boltz-2 primary validator, but
# with a different all-atom cofolder (Chai-1). Keeping designs that pass BOTH
# tools filters out model-specific artifacts — a standard consensus filter for
# de novo binders.
#
# Chai-1 is pip-installed (AVX2-safe, like Boltz) and runs on ESM embeddings
# without an MSA server (no compute-node egress needed): use_esm_embeddings=True,
# no --use-msa-server. Weights (~5 GB) download once to CHAI_DOWNLOADS_DIR.
#
# Required env vars:
#   MPNN_DIR      dir with seqs/*.fa (ProteinMPNN output)
#   TARGET_SEQ    target (receptor) amino-acid sequence
#   OUT_DIR       Chai output dir (e.g. .../04_chai)
#   CAP           max complexes to fold (0 = all); mirrors the Boltz cap so the
#                 two validators cover the same design set
#   CHAI_DOWNLOADS_DIR  persisted weights dir (exported by the caller)
set -eo pipefail

: "${MPNN_DIR:?set MPNN_DIR}"; : "${TARGET_SEQ:?set TARGET_SEQ}"
: "${OUT_DIR:?set OUT_DIR}"; : "${CAP:=0}"
: "${DIFFN_TIMESTEPS:=200}"; : "${RECYCLES:=3}"

mkdir -p "$OUT_DIR/fastas" "$OUT_DIR/preds"

# Build one Chai FASTA per designed sequence (target chain A + binder chain B),
# using the SAME design ids as the Boltz stage (<stem>_s<si>) so results join.
python - "$MPNN_DIR" "$OUT_DIR/fastas" "$TARGET_SEQ" "$CAP" <<'PYEOF'
import glob, os, sys
mpnn_dir, fa_dir, target_seq, cap = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
n = 0
paths = sorted(glob.glob(os.path.join(mpnn_dir, "seqs", "*.fa")))
for fa in paths:
    stem = os.path.basename(fa)[:-3]
    recs = open(fa).read().split(">")
    for si, rec in enumerate(recs[2:]):
        lines = rec.strip().split("\n")
        if len(lines) < 2:
            continue
        seq = "".join(lines[1:]).strip()
        if not seq:
            continue
        did = f"{stem}_s{si}".replace("/", "_").replace(".", "_")
        with open(os.path.join(fa_dir, f"{did}.fasta"), "w") as o:
            o.write(f">protein|name=target\n{target_seq}\n")
            o.write(f">protein|name=binder\n{seq}\n")
        n += 1
print(f"[chai] wrote {n} complex FASTAs")
PYEOF

# Apply the cap (keep first CAP by sorted design id) to match Boltz's design set.
if [ "$CAP" -gt 0 ]; then
  ls "$OUT_DIR"/fastas/*.fasta | sort | tail -n +$((CAP+1)) | xargs -r rm -f
  echo "[chai] capped to $CAP FASTAs"
fi

echo "[chai] folding complexes with Chai-1 (esm embeddings, no MSA)"
START=$SECONDS
python - "$OUT_DIR/fastas" "$OUT_DIR/preds" "$DIFFN_TIMESTEPS" "$RECYCLES" <<'PYEOF'
import glob, os, sys, json, numpy as np
from pathlib import Path
from chai_lab.chai1 import run_inference
fa_dir, pred_root, ndiff, nrec = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rows = []
for fa in sorted(glob.glob(os.path.join(fa_dir, "*.fasta"))):
    did = os.path.basename(fa)[:-6]
    odir = os.path.join(pred_root, did)
    if os.path.isdir(odir):  # run_inference refuses a non-empty dir
        import shutil; shutil.rmtree(odir)
    try:
        cands = run_inference(
            fasta_file=Path(fa), output_dir=Path(odir),
            num_trunk_recycles=nrec, num_diffn_timesteps=ndiff,
            seed=42, device="cuda:0", use_esm_embeddings=True,
        )
        # pick best sample by aggregate_score
        best = int(np.argmax([rd.aggregate_score.item() for rd in cands.ranking_data]))
        npz = os.path.join(odir, f"scores.model_idx_{best}.npz")
        d = np.load(npz)
        rows.append(dict(design=did,
                         chai_aggregate=float(d["aggregate_score"].reshape(-1)[0]),
                         chai_iptm=float(d["iptm"].reshape(-1)[0]),
                         chai_ptm=float(d["ptm"].reshape(-1)[0])))
    except Exception as e:
        rows.append(dict(design=did, chai_aggregate="", chai_iptm="", chai_ptm="", error=str(e)[:120]))
import csv
out = os.path.join(os.path.dirname(pred_root), "chai_scores.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["design","chai_aggregate","chai_iptm","chai_ptm","error"])
    w.writeheader()
    for r in rows: w.writerow(r)
print(f"[chai] scored {sum(1 for r in rows if r.get('chai_iptm') not in ('',None))}/{len(rows)} complexes -> {out}")
PYEOF
echo "[chai] runtime $((SECONDS-START))s"
NS=$(python -c "import csv,sys; rows=list(csv.DictReader(open('$OUT_DIR/chai_scores.csv'))); print(sum(1 for r in rows if r['chai_iptm'] not in ('',None)))" 2>/dev/null || echo 0)
echo "[chai] scored complexes: $NS"
[ "$NS" -gt 0 ] && echo "STAGE_CHAI_OK"
