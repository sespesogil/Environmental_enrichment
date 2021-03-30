#PBS -N EnhancER.EE
#PBS -l nodes=1:ppn=20
#PBS -l vmem=60g
#PBS -q parallel 
#PBS -l walltime=20:01:30
#PBS -joe /scripts/featureCounts

# Example of read counts retrieval over enhancers using featureCounts


module purge
module load python/2.7.9
module load samtools
module load subread
module load igvtools
module load bowtie/2.3.4.2



SAF=/GEP.combined.enhancers.mm10.saf
bams=/bamfiles
out=/enhancers

mkdir $out
cd $bams

find `pwd` -name \*.sort.bam -exec ls {} \; | sort >  $out/list.chip.fastq.txt
find `pwd` -name \*.sort.bam -printf '%f\n' | cut -d"." -f1 | sort | uniq > $out/filenames.chip.txt



list1=$out/list.chip.fastq.txt
list2=$out/filenames.chip.txt


paste -- "$list1" "$list2" |
while IFS=$'\t' read -r file1 file2 rest; do
mkdir $out/$file2
featureCounts -T 20 --ignoreDup -a $SAF -F SAF -O -o $out/$file2/$file2.enhancers.txt $file1
awk -F"\t" '{print $7}' $out/$file2/$file2.enhancers.txt > $out/$file2/$file2.enhancers.counts.txt
awk -F"\t" '{print $1}' $out/$file2/$file2.enhancers.txt > $out/features.enhancer.txt

  case $? in
    0) status='same';;
    1) status='different';;
    *) status='ERROR';;
  esac
  echo "$status $file1 $file2"
done

paste $out/*/*enhancers.counts.txt > $out/consensus.enhancers.counts.txt
paste $out/features.enhancer.txt $out/consensus.enhancers.counts.txt > $out/ENHANCERS.consensus.enhancers.counts.txt
