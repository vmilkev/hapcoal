#!/bin/sh
# Human chr22 validation: 1000 Genomes CEU LD decay vs SMC, SMC' and MaCS
# simulated with the HapMap II GRCh37 map and the Tennessen et al. (2012) CEU demography.
set -e
cd "$(dirname "$0")/.."
DATA=${DATA:-data}
OUT=${OUT:-out/chr22}
[ -s "$DATA/genetic_map_GRCh37_chr22.txt" ] && \
[ -s "$DATA/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" ] || sh scripts/fetch_chr22.sh "$DATA"
./bin/hapcoal chr22 --data "$DATA" --out "$OUT" "$@"
echo
column -t "$OUT/comparison.tsv" 2>/dev/null || cat "$OUT/comparison.tsv"
# PNG with matplotlib if installed, otherwise a dependency-free SVG
command -v python3 >/dev/null && python3 scripts/plot_ld.py "$OUT" || echo "(python3 not found: skipping plot)"
