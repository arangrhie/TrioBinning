#!/bin/bash

if [[ $# -lt 4 ]]; then
  echo "Usage: assign.sh hap1.hapmer.meryl hap2.hapmer.meryl input.fofn cutoff"
  echo "  hap1.hapmer.meryl  hap1 hapmer meryl dir path"
  echo "  hap2.hapmer.meryl  hap2 hapmer meryl dir path"
  echo "  input.fofn         input fastq files"
  echo "  cutoff             cutoff ratio for defining hap1 if hap1 / hap2 > cutoff"
  exit 0
fi

set -x
set -o pipefail

# Only one input supported for `-sequence`

hapmer1=$1
hapmer2=$2
input=$3
cutoff=$4

hap1=`basename $hapmer1 | sed 's/.hapmer.meryl//'`
hap2=`basename $hapmer2 | sed 's/.hapmer.meryl//'`

i=$SLURM_ARRAY_TASK_ID

if [[ -z $SLURM_ARRAY_TASK_ID ]]; then
  echo "No SLURM_ARRAY_TASK_ID variable detected."
  echo "Default to run on the 1st line of input.fofn"
  i=1
fi

# Only one input supported for `-sequence`
fq=`sed -n ${i}p $input`
name=`basename $fq`
echo "Running on ${i}th line of $input: $fq"

# $4/$3 is the normalized count for hap1, $6/$5 is the normalized count for hap2
# score 1 is always larger or equal to score 2
# if score 2 is 0, score 1 is larger than 0, assign to hap1
# if score 2 is larger than 0, and score 1/score 2 > cutoff, assign to hap
# otherwise, assign to "unknown"

input=`echo $input | sed 's/.fofn//g'`
out=$input.assigned.$i.$cutoff

if [[ -s $out ]]; then
  echo "$out will be kept as $out.bk"
  mv $out $out.bk
fi
meryl-lookup -existence -sequence $fq -mers $hapmer1 $hapmer2 |\
  awk -v hap1=$hap1 -v hap2=$hap2 -v cutoff=$cutoff \
    '{sco1=$4/$3; sco2=$6/$5; hap=hap1; \
      if(sco2>sco1) {tmp=sco1; sco1=sco2; sco2=tmp; hap=hap2;}; \
      if((sco2==0 && sco1>0) || (sco2>0 && (sco1/sco2>cutoff))) { print $0"\t"hap } \
      else {print $0"\tunk"} \
    }' >> $out

#module load python
#python $tools/TrioBinning/bin_reads.py $fq $out $cutoff

set +x
