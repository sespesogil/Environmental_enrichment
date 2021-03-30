#!/bin/bash
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -W 30:00
#BSUB -J diffReps_ATACseq_cortex
#BSUB -o diffReps_XOvWTmg.o
#BSUB -e diffReps_XOvWTmg.e

# Example of DiffReps run

OUT=/analysis/diffReps/cortex/ATACseq

cd $OUT

module load diffreps/1.55.4 R bedtools

glia=/bamfiles/NeuN-/H3K79me2
neurons=/bamfiles/NeuNpos/H3K79me2
cortex=/bamfiles/cortex/ATACseq

bedtools bamtobed -i $cortex/WEN1.dedup_removed.sort.bam > $cortex/WEN1.dedup_removed.bed
bedtools bamtobed -i $cortex/WEN2.dedup_removed.sort.bam > $cortex/WEN2.dedup_removed.bed
bedtools bamtobed -i $cortex/WNN1.dedup_removed.sort.bam > $cortex/WNN1.dedup_removed.bed
bedtools bamtobed -i $cortex/WNN2.dedup_removed.sort.bam > $cortex/WNN2.dedup_removed.bed

bedtools sort -i $cortex/WEN1.dedup_removed.bed
bedtools sort -i $cortex/WEN2.dedup_removed.bed
bedtools sort -i $cortex/WNN1.dedup_removed.bed
bedtools sort -i $cortex/WNN2.dedup_removed.bed

diffReps.pl -tr $cortex/WEN1.dedup_removed.bed $cortex/WEN2.dedup_removed.bed -co $cortex/WNN1.dedup_removed.bed $cortex/WNN2.dedup_removed.bed -re diffReps_ATACseq_EE_vs_CTL_1kb -meth gt -chrlen mm10.chrom.sizes -pval 0.001 -frag 150 -window 1000
 


