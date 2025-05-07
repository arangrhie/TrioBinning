#!/bin/bash

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 mat.hapmer.meryl pat.hapmer.meryl input.fofn [cutoff]"
  echo ""
  echo "Assign each read to hap1, hap2 or Unassigned."
  echo "Required:"
  echo "  mat.hapmer.meryl"
  echo "  pat.hapmer.meryl"
  echo "  input.fofn"
  echo "  cutoff            DEFAULT=4.0"
  exit 0
fi

mat=$1
pat=$2
input=$3
cutoff=$4

if [[ -z $cutoff ]]; then
  cutoff=4.0
fi

mkdir -p logs

cpus=12
mem=24g
name=assign
script=$tools/TrioBinning/assign.sh
args="$mat $pat $input $cutoff"
partition=norm
walltime=24:00:00
path=`pwd`
extra="--array=1-$(wc -l $input | cut -f1 -d' ')"
log=logs/$name.%A_%a.log

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args

