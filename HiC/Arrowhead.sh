#!/bin/bash
#BSUB -n 6
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=15000]"
#BSUB -W 100:00
#BSUB -J juicer_tools_arrowhead
#BSUB -o 70i79854di66267806965718.o
#BSUB -e 5602465i670d946i5677418.e


# TAD calling for Chrom3D pipeline 

module purge
module load juicer/1.5_cpu


jarcuda=/juicer/1.5/cpu/scripts/juicer_tools.1.7.6_jcuda.0.8.jar

hic=EE.allValidPairs.hic
dir=Arrowhead/EE
mkdir $dir

java -jar $jarcuda arrowhead -r 10000 $hic $dir/contact_domains_list_10kb --ignore_sparsity
java -jar $jarcuda arrowhead -r 10000 $hic $dir/contact_domains_list_10kb_VC -k VC --ignore_sparsity
java -jar $jarcuda arrowhead -r 10000 $hic $dir/contact_domains_list_10kb_VC_SQRT -k VC_SQRT --ignore_sparsity
java -jar $jarcuda arrowhead -r 10000 $hic $dir/contact_domains_list_10kb_KR -k KR --ignore_sparsity

java -jar $jarcuda arrowhead -r 25000 $hic $dir/contact_domains_list_25kb --ignore_sparsity
java -jar $jarcuda arrowhead -r 25000 $hic $dir/contact_domains_list_25kb_VC -k VC --ignore_sparsity
java -jar $jarcuda arrowhead -r 25000 $hic $dir/contact_domains_list_25kb_VC_SQRT -k VC_SQRT --ignore_sparsity
java -jar $jarcuda arrowhead -r 25000 $hic $dir/contact_domains_list_25kb_KR -k KR --ignore_sparsity

java -jar $jarcuda arrowhead -r 50000 $hic $dir/contact_domains_list_50kb --ignore_sparsity
java -jar $jarcuda arrowhead -r 50000 $hic $dir/contact_domains_list_50kb_VC -k VC --ignore_sparsity
java -jar $jarcuda arrowhead -r 50000 $hic $dir/contact_domains_list_50kb_VC_SQRT -k VC_SQRT --ignore_sparsity
java -jar $jarcuda arrowhead -r 50000 $hic $dir/contact_domains_list_50kb_KR -k KR --ignore_sparsity
