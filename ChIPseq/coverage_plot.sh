#!/bin/bash
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=18000]"
#BSUB -W 60:00
#BSUB -J deepTools_v4
#BSUB -o o_deeptools_0204183.o
#BSUB -e e_deeptools_02049183.e

# Example of coverage plots using DEEPTOOLS 

module load python/2.7.14
module load py_packages


DEEPTOOLS=/programs/deepTools2.0/bin
OUT=/deepTools/02_09_18_ATACseq_cov_cortex_differential_interactions_GAINS

BAM=/bamfiles/cortex

$DEEPTOOLS/bamCoverage -b $BAM/WEN1.dedup_removed.sort.bam -p 4 -o $OUT/WEN1.dedup_removed.bw  --normalizeUsingRPKM
$DEEPTOOLS/bamCoverage -b $BAM/WEN2.dedup_removed.sort.bam -p 4 -o $OUT/WEN2.dedup_removed.bw  --normalizeUsingRPKM
$DEEPTOOLS/bamCoverage -b $BAM/WNN1.dedup_removed.sort.bam -p 4 -o $OUT/WNN1.dedup_removed.bw  --normalizeUsingRPKM
$DEEPTOOLS/bamCoverage -b $BAM/WNN2.dedup_removed.sort.bam -p 4 -o $OUT/WNN2.dedup_removed.bw  --normalizeUsingRPKM


BED=/analysis/bedfiles

$DEEPTOOLS/computeMatrix reference-point -S $OUT/WEN1.dedup_removed.bw $OUT/WEN2.dedup_removed.bw $OUT/WNN1.dedup_removed.bw $OUT/WNN2.dedup_removed.bw -R $BED/sig_anchor_gains.bed  -a 1000 -b 1000 -o $OUT/matrix_reference.gz
$DEEPTOOLS/computeMatrix scale-regions -S $OUT/WEN1.dedup_removed.bw $OUT/WEN2.dedup_removed.bw $OUT/WNN1.dedup_removed.bw $OUT/WNN2.dedup_removed.bw -R $BED/sig_anchor_gains.bed -a 1000 -b 1000 -o $OUT/matrix_scaled.gz

$DEEPTOOLS/plotProfile -m $OUT/matrix_reference.gz -out $OUT/matrix_reference_profile_anchors.pdf --perGroup
$DEEPTOOLS/plotHeatmap -m $OUT/matrix_reference.gz -out $OUT/matrix_reference_anchors.pdf

$DEEPTOOLS/plotProfile -m $OUT/matrix_scaled.gz -out $OUT/matrix_scaled_profile_anchors.pdf --perGroup
$DEEPTOOLS/plotHeatmap -m $OUT/matrix_scaled.gz -out $OUT/matrix_scaled_anchors.pdf





