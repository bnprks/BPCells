set -euo pipefail

INPUT_FILES="/oak/stanford/groups/wjg/bparks/BPCells/01_raw_data/ENCODE/samples/snyder_*_PANC/fragments.tsv.gz"

GENOME="hg38"
MIN_FRAGS="2000"
MIN_TSS="5"

OUTPUT_DIR="/oak/stanford/groups/wjg/bparks/BPCells/04_data/test_real_data/bpcells_snyder_panc"
PEAK_SET="/oak/stanford/groups/wjg/bparks/BPCells/01_raw_data/peaksets/ENCFF305BGH.shuffle.bed"

THREADS="8"

SINGULARITY_COMMAND="singularity run --no-home --cleanenv --env OMP_NUM_THREADS=1  $OAK/bparks/BPCells/02_config/envs/archr_bpcells.sif"
$SINGULARITY_COMMAND Rscript bpcells_archr_compatible.R "$INPUT_FILES" "$GENOME" "$MIN_FRAGS" "$MIN_TSS" "$OUTPUT_DIR" "$PEAK_SET" "$THREADS"