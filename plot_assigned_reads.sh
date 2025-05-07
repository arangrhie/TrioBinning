#!/bin/sh

if [[ $# -lt 3 ]]; then
  echo "Usage: plot_assigned_reads.sh input.fofn hap1 hap2 [cutoff]"
  echo "  hap1    hap1 name on the plot e.g. Mat"
  echo "  hap2    hap2 name on the plot e.g. Pat"
  echo "  out     output file name prefix"
  echo "  cutoff  cutoff ratio for defining hap1 if norm(hap1) / norm(hap2) > cutoff. DEFAULT=4.0"
  exit 0
fi

input=$1
hap1=$2
hap2=$3
cutoff=$4
if [[ -z $cutoff ]]; then
  cutoff=4.0
fi
min=100

LEN=`wc -l $input | cut -f1 -d' '`
echo $LEN

input=`echo $input | sed 's/.fofn//g'`
out=$input.assigned.$cutoff

if [[ -s $out ]]; then
  echo "$out will be kept as $out.bk"
  mv $out $out.bk
fi

for i in $(seq 1 $LEN); do
  in=$input.assigned.$i.$cutoff
  cat $in |
    awk -v hap1=$hap1 -v hap2=$hap2 -v cutoff=$cutoff '{sco1=$4/$3; sco2=$6/$5; hap=hap1; \
      if(sco2>sco1) {tmp=sco1; sco1=sco2; sco2=tmp; hap=hap2;}; \
      if((sco2==0 && sco1>0) || (sco2>0 && (sco1/sco2>cutoff))) { print $0"\t"hap } \
      else if((sco2>0 && sco1/sco2<=cutoff) && $4>100 && $6>100) { print $0"\tmix" } \
      else {print $0"\tunk"} \
     }' >> $out
done

echo -e "Hap\tReadId\tMat\tPat\tReadLen" > $out.$min.blob
awk -v min=$min '$4>min && $6>min {print $NF"\t"$1"\t"$4"\t"$6"\t"$2}' $out >> $out.$min.blob

module load R
Rscript $tools/merqury/plot/plot_blob.R -f $out.$min.blob -o $out.$min.blob
