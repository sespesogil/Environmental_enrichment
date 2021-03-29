#!/bin/bash
#BSUB -q premium
#BSUB -n 3
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=15000]"
#BSUB -W 100:00
#BSUB -J juicer_tools_CTL2
#BSUB -o 670i7854i6626780465676918.o
#BSUB -e 5602465670946i56766718.e

# example to compute eigen verctor values for chr19 
module purge
module load juicer/1.5_cpu
module unload glibc/2.14


java -jar /scripts/juicer_tools.1.7.6_jcuda.0.8.jar eigenvector NONE CTL_rep2_allValidPairs.hic 19 BP 1000000 eigen_chr19_NONE_1000000.txt

java -jar /scripts/juicer_tools.1.7.6_jcuda.0.8.jar eigenvector VC CTL_rep2_allValidPairs.hic 19 BP 1000000 eigen_chr19_VC_1000000.txt

java -jar /scripts/juicer_tools.1.7.6_jcuda.0.8.jar eigenvector VC_SQRT CTL_rep2_allValidPairs.hic 19 BP 1000000 eigen_chr19_VC_SQRT_1000000.txt

java -jar /scripts/juicer_tools.1.7.6_jcuda.0.8.jar eigenvector KR /CTL_rep2_allValidPairs.hic 19 BP 1000000 eigen_chr19_KR_1000000.txt

