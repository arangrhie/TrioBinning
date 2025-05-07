#!/bin/bash

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 input.fofn cutoff"
  echo ""
  echo "Bin each read to hap1, hap2 or Unassigned."
  echo "Required:"
  echo "  input.fofn"
  echo "  cutoff            DEFAULT=4.0"
  exit 0
fi

input=$1
cutoff=$2

if [[ -z $cutoff ]]; then
  cutoff=4.0
fi

mkdir -p logs

cpus=4
mem=4g
name=bin
script=$tools/TrioBinning/bin_reads.sh
args="$input"
partition=norm
walltime=24:00:00
path=`pwd`
extra="--array=1-$(wc -l $input | cut -f1 -d' ') --dependency=afterok:$(cat assign.jid | awk '{print $NF}')"
log=logs/$name.%A_%a.log

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > bin.jid
set +x


