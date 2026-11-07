#!/bin/sh
# Output regression test: runs a fixed set of hapcoal commands and compares the
# SHA-256 checksums of all (deterministic) outputs with tests/golden.sha256.
#
#   sh tests/regression.sh check     compare with the recorded checksums (default)
#   sh tests/regression.sh record    re-record the checksums (only after a deliberate
#                                    change of results!)
#
# Requires the chr22 data (make data). Outputs go to out/regression.
#
# The reference checksums were recorded with Apple clang 21 (arm64). Other
# compilers may contract floating-point multiply-adds differently (e.g. GCC's
# default -ffp-contract=fast), which can change the last bit of a SNP position
# and thereby move a single SNP pair to the neighbouring LD bin; haplotypes,
# maps and fitted demographies are unaffected. Compare with the same compiler.
set -e
cd "$(dirname "$0")/.."
MODE=${1:-check}
BIN=${HAPCOAL:-./bin/hapcoal}
O=out/regression
DATA=data
rm -rf "$O"
mkdir -p "$O"

run() { "$BIN" "$@" --threads 1 --quiet; }

# 1. SMC, two chromosomes, hotspot maps, text output, LD matrix, LD-based map
run simulate --model smc --n-ind 50 --chr 1:3000000 --chr 2:2000000 --hotspots --text \
    --estimate-map --ld-matrix --ld-matrix-dist 20000 --seed 11 --out "$O/r1"
# 2. SMC', ploidy 3, real map region, bottleneck demography
run simulate --model smcp --n-ind 40 --ploidy 3 --map "$DATA/genetic_map_GRCh37_chr22.txt" \
    --chr 22:20000001-24000000 --epoch 0:20000 --epoch 300:2000 --epoch 1500:12000 --seed 12 --out "$O/r2"
# 3. MaCS, haplotype count, demography, genotype LD, dense LD region
run simulate --model macs --n-hap 100 --chr 1:5000000 --epoch 0:40000 --epoch 200:3000 --epoch 2500:12000 \
    --history 50000 --ld-genotype --ld-region 1:1000000-1100000 --seed 13 --out "$O/r3"
# 4. MaCS with a long history on the real map
run simulate --model macs --n-ind 30 --map "$DATA/genetic_map_GRCh37_chr22.txt" --chr 22:30000001-33000000 \
    --history 1000000 --seed 14 --out "$O/r4"
# 5. observed LD from the 1000 Genomes VCF
"$BIN" ld --threads 1 --vcf "$DATA/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
    --panel "$DATA/integrated_call_samples_v3.20130502.ALL.panel" --pop CEU --region 16050000-20000000 \
    --out "$O/obs" 2>/dev/null
# 6. chr22 test case with a short demography fit (SMC' only)
"$BIN" chr22 --data "$DATA" --out "$O/chr22" --models smcp --fit-evals 12 --threads 1 >/dev/null 2>&1

# Checksums of all outputs except logs, plots and timing columns.
SUM=$(command -v sha256sum >/dev/null && echo "sha256sum" || echo "shasum -a 256")
cut -f1-7,9- "$O/chr22/comparison.tsv" > "$O/chr22/comparison.notime.tsv"
( cd "$O" && find . -type f ! -name '*.log' ! -name '*.svg' ! -name '*.png' ! -name 'comparison.tsv' ! -name 'current.sha256' \
    | LC_ALL=C sort | xargs $SUM ) > "$O/current.sha256"

if [ "$MODE" = record ]; then
    cp "$O/current.sha256" tests/golden.sha256
    echo "recorded $(wc -l < tests/golden.sha256) checksums in tests/golden.sha256"
else
    if diff tests/golden.sha256 "$O/current.sha256" >/dev/null; then
        echo "regression: all $(wc -l < tests/golden.sha256) outputs identical to the recorded results"
    else
        echo "regression: OUTPUTS DIFFER from the recorded results:"
        diff tests/golden.sha256 "$O/current.sha256" || true
        exit 1
    fi
fi
