#!/usr/bin/env bash
# Clean-room reproduction of Table 2 from public downloads (see README_pipeline.md, "Reproduce Table 2 from scratch").
#   bash reproduce_table2.sh            # run in a copy of this folder; needs ~25 GB free disk, ~8 GB RAM, internet
# Steps are check-pointed (data/derived/.done_<step>); re-running resumes after the last completed step.
# Manual prerequisite: the three DepMap-portal files listed in data_manifest.csv (retrieval = manual) saved under data/.
set -euo pipefail
cd "$(dirname "$0")"
export GYN_FIDELITY_PATHS="$PWD/config_paths.cleanroom.json"
LOG="$PWD/reproduce_table2.log"; mkdir -p data/derived
step(){ local name=$1; shift
  if [ -e "data/derived/.done_$name" ]; then echo "$(date '+%F %T') SKIP $name (done)" | tee -a "$LOG"; return; fi
  local t0=$(date +%s); echo "$(date '+%F %T') START $name: $*" | tee -a "$LOG"
  "$@" 2>&1 | tee -a "$LOG"
  echo "$(date '+%F %T') END $name ($(( $(date +%s)-t0 )) s)" | tee -a "$LOG"; touch "data/derived/.done_$name"; }
step fetch      python3 fetch_public_data.py --data data
step hcmi_rna   python3 00e_build_hcmi_rna_and_gene_annotation.py
step purity     python3 00f_build_purity_filtered_genes.py
step reference  python3 01_build_reference.py
step models     python3 02_build_model_matrix.py
step run_26Q3   python3 run_all.py --depmap-expr data/depmap26q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv --tag 26Q3 --expect-models 2130 --no-outputs
step compare    python3 tests/compare_table2.py --results results/26Q3 --out reproduce_table2_check.json
echo "$(date '+%F %T') Table 2 reproduced cell for cell (reproduce_table2_check.json)" | tee -a "$LOG"
