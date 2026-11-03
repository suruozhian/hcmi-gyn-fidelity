#!/usr/bin/env bash
# Clean-room reproduction of every analysis reported in the manuscript from public data (see README_pipeline.md,
# "Reproducibility scope"). Run inside supp/code of a copy of the manuscript package (manuscript.md two levels up):
#   bash reproduce_all.sh [published_package_root]     # default ../.. ; needs ~27 GB disk, internet
# Steps are check-pointed (data/derived/.done_<step>); a re-run resumes after the last completed step.
# Manual prerequisite: the DepMap-portal files of data_manifest.csv (retrieval = manual) saved under data/.
set -euo pipefail
cd "$(dirname "$0")"
PUB="${1:-../..}"
export GYN_FIDELITY_PATHS="$PWD/config_paths.cleanroom.json"
LOG="$PWD/reproduce_all.log"; mkdir -p data/derived
step(){ local name=$1; shift
  if [ -e "data/derived/.done_$name" ]; then echo "$(date '+%F %T') SKIP $name (done)" | tee -a "$LOG"; return; fi
  local t0=$(date +%s); echo "$(date '+%F %T') START $name: $*" | tee -a "$LOG"
  "$@" 2>&1 | tee -a "$LOG"
  echo "$(date '+%F %T') END $name ($(( $(date +%s)-t0 )) s)" | tee -a "$LOG"; touch "data/derived/.done_$name"; }
D26=data/depmap26q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
# inputs
step fetch           python3 fetch_public_data.py --data data
step hcmi_rna        python3 00e_build_hcmi_rna_and_gene_annotation.py
step hcmi_dna        python3 00h_build_hcmi_dna.py
step purity          python3 00f_build_purity_filtered_genes.py
step reference       python3 01_build_reference.py
step reference_sym   python3 01_build_reference.py --ids symbol
step models          python3 02_build_model_matrix.py
step models_sym      python3 00g_build_symbol_model_matrix.py
# molecular derivations (Fig. 3, Fig. 4, Table 3, Tables S9/S10, DNA columns of Table S2, Supplementary Fig. S7)
step molecular_hcmi  python3 04_derive_hcmi_molecular.py
step endometrial     python3 05_derive_endometrial_classes.py
step fig_surrogate   python3 06_fig_hcmi_surrogate.py
# analyses: headline (DepMap 26Q3, Entrez) and sensitivity runs (DepMap 24Q4; gene symbols); each run includes the
# reference-composition (non-gyn 40/project; gyn ~200/project) and equal-depth (40/project x 25) analyses
step run_26Q3        python3 run_all.py --depmap-expr $D26 --tag 26Q3 --expect-models 2130 --no-outputs
step run_24Q4        python3 run_all.py --tag 24Q4 --no-outputs
step run_26Q3_symbol python3 run_all.py --depmap-expr $D26 --tag 26Q3_symbol --gene-id symbol --expect-models 2130 --no-outputs
step run_24Q4_symbol python3 run_all.py --tag 24Q4_symbol --gene-id symbol --no-outputs
step identity_26Q3   python3 identity_test.py --tag 26Q3 --depmap-expr $D26
step identity_24Q4   python3 identity_test.py --tag 24Q4
step identity_26Q3s  python3 identity_test.py --tag 26Q3_symbol --gene-id symbol --depmap-expr $D26
step histotype_ref   bash -c "python3 03_histotype_reference_checks.py --out results/26Q3/histotype_reference_checks.json && cp results/26Q3/histotype_reference_checks.json results/24Q4/ && python3 03_histotype_reference_checks.py --gene-id symbol --out results/26Q3_symbol/histotype_reference_checks.json && cp results/26Q3_symbol/histotype_reference_checks.json results/24Q4_symbol/"
step regression_syn  python3 tests/test_gene_alignment.py --synthetic --out results/26Q3/regression_test_synthetic.json
# the real-data regression test FAILS by design of its pre-specified criterion C2 (published result; README "Gene alignment and
# regression test"): its exit code 1 is recorded, not treated as an error; its output is compared with the published file
step regression_real bash -c "python3 tests/test_gene_alignment.py --real --tag 26Q3 --depmap-expr $D26 --out results/26Q3/regression_test_real.json; echo \"exit code \$? (expected 1: pre-specified criterion C2 not met, as published)\""
step regression_diag python3 tests/diagnose_regression.py --tag 26Q3 --depmap-expr $D26 --cases 5:0,20:0,36:4 --out results/26Q3/regression_test_real_diagnostic.json
# outputs: figures, tables, source data and number audit written into pkg_out/ (a package copy), then compared
step outputs         bash -c "mkdir -p pkg_out && cp '$PUB/manuscript.md' pkg_out/manuscript.md && python3 make_outputs.py --tag 24Q4 && python3 make_outputs.py --tag 26Q3_symbol && python3 make_outputs.py --tag 26Q3 --publish"
# comparisons: exit codes are recorded; the summary below fails if any comparison found a difference
step compare_table2  python3 tests/compare_table2.py --results results/26Q3 --out reproduce_table2_check.json
step compare_molecular bash -c "python3 tests/compare_molecular.py --out molecular_rederivation_check.json; echo \"exit code \$?\""
step compare_all     bash -c "python3 tests/compare_all.py --fresh pkg_out --published '$PUB' --out reproduce_all_check.json; echo \"exit code \$?\""
python3 - <<'PY' | tee -a "$LOG"
import json, sys
m=json.load(open('molecular_rederivation_check.json'))['summary']['all_match']; c=json.load(open('reproduce_all_check.json'))['identical']
print('molecular intermediates match the shipped copies (tests/compare_molecular.py):',m,'; outputs identical to the published package:',c)
sys.exit(0 if m and c else 1)
PY
