#!/bin/bash

if [[ $# -lt 2 ]]; then
  echo "Usage: bin_reads.sh input.fofn cutoff"
  echo "  input.fofn  input fastq files"
  exit 0
fi

set -x
set -o pipefail

input=$1
cutoff=$2

i=$SLURM_ARRAY_TASK_ID

if [[ -z $SLURM_ARRAY_TASK_ID ]]; then
  echo "No SLURM_ARRAY_TASK_ID variable detected."
  echo "Default to run on the 1st line of input.fofn"
  i=1
fi

# Only one input supported for `-sequence`
fq=`sed -n ${i}p $input`
name=`basename $fq`

input=`echo $input | sed 's/.fofn//g'`
out=$input.assigned.$i.$cutoff

echo "Running on ${i}th line of $input: $name"

module load python
python $tools/TrioBinning/bin_reads.py $fq $out $cutoff
