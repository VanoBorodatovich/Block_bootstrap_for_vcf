#!/usr/bin/bash
set -e

### INPUT
VCF="$1"
N_BOOT=100
THREADS=10
REC_RATE=1e-6
PROJECTION=14,6
# put projection in the same order as your popmap
POPMAP=pop_norw_2pop
OUTDIR="boot_data"

### PART 1 - GET SCAFFOLDS ACTUALLY PRESENT IN VCF
echo -e "chr\tlength" >> scaff_lengths_filtered.txt
bcftools view -h "$VCF" | \
grep '^##contig' | \
sed -E 's/.*ID=([^,]+),length=([0-9]+).*/\1\t\2/' > scaff_lengths.txt
 get list of scaffolds present in VCF
bcftools query -f '%CHROM\n' "$VCF" | sort -u > keep.txt
 filter your length table (chr\tlength)
grep -w -F -f keep.txt scaff_lengths.txt >> scaff_lengths_filtered.txt
rm keep.txt scaff_lengths.txt


### PART 2 - CREATE BOOTSTRAPS COORDINATES
Rscript - <<EOF
library(data.table)

calc_distance <- function(r, p=0.5) {
  -log(1 - p) / r
}

split_scaffold_equal <- function(L, d) {
  k <- max(1, round(L / d))
  base_size <- L %/% k
  remainder <- L %% k
  sizes <- rep(base_size, k)
  if (remainder > 0) sizes[1:remainder] <- sizes[1:remainder] + 1
  starts <- cumsum(c(0, sizes[-length(sizes)]))
  ends <- cumsum(sizes)
  data.table(start = starts, end = ends)
}

make_blocks <- function(scaffold_lengths, target_distance) {
  blocks <- scaffold_lengths[,
    split_scaffold_equal(length, target_distance),
    by = chr]
  blocks[, block_id := .I]
  setcolorder(blocks, c("block_id", "chr", "start", "end"))
  blocks[]
}

make_bootstrap_coords <- function(blocks_dt, n_boot, out_prefix="boot") {
  n_blocks <- nrow(blocks_dt)
  for (i in seq_len(n_boot)) {
    sampled <- blocks_dt[sample(.N, n_blocks, replace = TRUE)]
    sampled[, coord := paste0(chr, ":", start, "-", end)]
    fwrite(sampled[, .(coord)],
           file = sprintf("%s_%03d.txt", out_prefix, i),
           col.names = FALSE)
  }
}

# ---- pipeline ----
chr_l <- fread("scaff_lengths_filtered.txt")

dist <- calc_distance($REC_RATE, 0.5)
cat("distance for block bootstrap is ", dist, "\n")
blocks_dt <- make_blocks(chr_l, dist)

make_bootstrap_coords(blocks_dt, $N_BOOT, "boot")
EOF

### PART 3 - CREATE BOOTSRAP VCFS PARALLELY
ls boot_*.txt | parallel --jobs "$THREADS" 'bcftools view -h "$VCF" > {.}.vcf'
# define function, that would create vcf with duplicates positions
process_boot() {
  local file="$1"
  local out="${file%.txt}.vcf"
  for region in $(cat "$file"); do
    bcftools view -H -r "$region" "$VCF" >> "$out"
  done
}
export -f process_boot
export VCF
parallel --jobs $THREADS --verbose process_boot ::: boot_*.txt

### PART 4 - RUN EASYSFS ON EACH BOOTSTRAP PARALLELY
ls boot_*vcf | parallel --jobs "$THREADS" --verbose "easySFS.py -a -i {} --proj $PROJECTION -o {.}.sfs -p $POPMAP"
mkdir -p $OUTDIR
for f in boot*sfs; do mv ${f}/dadi/kola-norw.sfs ${OUTDIR}/${f}; done
