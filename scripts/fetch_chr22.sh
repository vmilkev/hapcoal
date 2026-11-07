#!/bin/sh
# Download the human chr22 test data (GRCh37):
#   * 1000 Genomes phase 3 chr22 phased VCF (~206 MB) + sample panel
#   * HapMap II GRCh37 recombination map (chr22 extracted from the archive)
set -e
DATA=${1:-data}
mkdir -p "$DATA"
cd "$DATA"

KG=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
VCF=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
PANEL=integrated_call_samples_v3.20130502.ALL.panel
HAPMAP=https://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz

[ -s "$PANEL" ] || curl -fL -o "$PANEL" "$KG/$PANEL"
[ -s "$VCF" ]   || curl -fL -o "$VCF.part" "$KG/$VCF" && mv -f "$VCF.part" "$VCF" 2>/dev/null || true
if [ ! -s genetic_map_GRCh37_chr22.txt ]; then
    curl -fL -o hapmap.tar.gz "$HAPMAP"
    tar -xzf hapmap.tar.gz genetic_map_GRCh37_chr22.txt
    rm -f hapmap.tar.gz
fi
ls -l "$PANEL" "$VCF" genetic_map_GRCh37_chr22.txt
