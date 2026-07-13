#!/bin/bash
#SBATCH --job-name=nkg7_dscalib
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=mlgpu_short
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.out
#
# ds_calib.sbatch — measure Boltz-2 ds=1 vs ds=5 scoring variance on 24 designs
# spanning the production composite-score range. Re-folds each design's EXISTING
# production YAML at ds=1 and ds=5, then compares scores + rank shift.
set -eo pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate boltz

WS=/lustre/scratch/data/dcorvino_hpc-nkg7_binder_design
CAL=$WS/runs/ds_calib
mkdir -p $CAL/yamls $CAL/ds1 $CAL/ds5

# ID list: arm:design (filled by render)
DESIGNS="human_ecl1:human_ecl1_472_s8 human_ecl1:human_ecl1_326_s11 human_ecl1:human_ecl1_140_s15 human_ecl1:human_ecl1_471_s12 human_ecl1:human_ecl1_444_s9 human_ecl1:human_ecl1_242_s7 human_ecl2:human_ecl2_482_s9 human_ecl2:human_ecl2_59_s6 human_ecl2:human_ecl2_228_s11 human_ecl2:human_ecl2_389_s0 human_ecl2:human_ecl2_27_s0 human_ecl2:human_ecl2_30_s6 mouse_ecl1:mouse_ecl1_95_s8 mouse_ecl1:mouse_ecl1_388_s6 mouse_ecl1:mouse_ecl1_249_s8 mouse_ecl1:mouse_ecl1_486_s1 mouse_ecl1:mouse_ecl1_19_s2 mouse_ecl1:mouse_ecl1_199_s7 mouse_ecl2:mouse_ecl2_196_s5 mouse_ecl2:mouse_ecl2_251_s2 mouse_ecl2:mouse_ecl2_447_s5 mouse_ecl2:mouse_ecl2_144_s0 mouse_ecl2:mouse_ecl2_114_s9 mouse_ecl2:mouse_ecl2_258_s15"

# gather the 24 production YAMLs into one dir
for row in $DESIGNS; do
  arm=${row%%:*}; d=${row##*:}
  cp $WS/runs/production/$arm/03_boltz/yamls/$d.yaml $CAL/yamls/$d.yaml
done
echo "[dscalib] gathered $(ls $CAL/yamls/*.yaml | wc -l) YAMLs"

# fold all at ds=1 then ds=5 (same inputs, different diffusion_samples)
echo "===== ds=1 ====="
boltz predict $CAL/yamls --out_dir $CAL/ds1 --diffusion_samples 1 \
  --recycling_steps 3 --output_format pdb --no_kernels --override
echo "===== ds=5 ====="
boltz predict $CAL/yamls --out_dir $CAL/ds5 --diffusion_samples 5 \
  --recycling_steps 3 --output_format pdb --no_kernels --override

# compare: pull ipTM/complex_plddt from each confidence json, join on design
python - "$CAL" <<'PYEOF'
import json, glob, os, sys, csv
cal = sys.argv[1]
def scores(root):
    out={}
    for f in glob.glob(os.path.join(root,"**","confidence_*_model_0.json"), recursive=True):
        d=json.load(open(f))
        name=os.path.basename(f).replace("confidence_","").replace("_model_0.json","")
        out[name]=dict(iptm=d.get("iptm"), plddt=d.get("complex_plddt"),
                       ptm=d.get("ptm"), conf=d.get("confidence_score"))
    return out
s1=scores(os.path.join(cal,"ds1","boltz_results_yamls","predictions"))
s5=scores(os.path.join(cal,"ds5","boltz_results_yamls","predictions"))
# fallback glob if structure differs
if not s1: s1=scores(os.path.join(cal,"ds1"))
if not s5: s5=scores(os.path.join(cal,"ds5"))
rows=[]
for name in sorted(set(s1)|set(s5)):
    a=s1.get(name,{}); b=s5.get(name,{})
    rows.append(dict(design=name, iptm_ds1=a.get("iptm"), iptm_ds5=b.get("iptm"),
                     plddt_ds1=a.get("plddt"), plddt_ds5=b.get("plddt"),
                     conf_ds1=a.get("conf"), conf_ds5=b.get("conf")))
with open(os.path.join(cal,"ds_compare.csv"),"w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"[dscalib] wrote ds_compare.csv with {len(rows)} designs")
PYEOF
echo "DSCALIB_DONE"
