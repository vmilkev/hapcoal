#!/bin/sh
# Whole-genome benchmark: 5000 diploid individuals (10,000 haplotypes) on the 22
# GRCh37 autosomes. Output goes to $OUT (default out/bench, removed by make rclean).
# The ~18 GB haplotype file is deleted after the run
# unless KEEP=1.
set -e
cd "$(dirname "$0")/.."
OUT=${OUT:-out/bench}
MODEL=${MODEL:-macs}
mkdir -p "$OUT"
TIME="/usr/bin/time"
case "$(uname -s)" in
    Darwin) TIME="/usr/bin/time -l" ;;
    Linux)  TIME="/usr/bin/time -v" ;;
esac
$TIME ./bin/hapcoal simulate --model "$MODEL" --n-ind 5000 --ploidy 2 --Ne 10000 --mu 1.25e-8 \
    --chromosomes data/grch37_autosomes.tsv --out "$OUT/genome_$MODEL" "$@"
ls -l "$OUT/genome_$MODEL.hap.bin"
[ "${KEEP:-0}" = 1 ] || rm -f "$OUT/genome_$MODEL.hap.bin"
