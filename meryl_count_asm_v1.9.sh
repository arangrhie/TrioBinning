#!/bin/bash

set -eu

me=$(basename "$0")
########################################################################################################################
# helper functions
show_help() {
cat << EOF
Takes in meryl databases and diploid assemblies, and outputs counts file for making blob plots.
Assumes ploidy == 2.

Syntax
  ${me} \\
      [Options] \\
      <mer_size> \\
      <hapA.meryl.db> <hapB.meryl.db> \\
      <asmA.fasta> <asmB.fasta>

  Options:
    -h or --help:
      display help and exit
    -p or --prefix:
      prefix to output counts files

  Mandatory Positional Arguments:
    mer_size:           selected k for counting
    hap(A|B).meryl.db:  path to meryl database directory for haplotype (A|B)
    asm(A|B).fasta:     path to assembly fasta for haplotype (A|B)
EOF
}

throw_error() {
    echo "$1" >&2
    exit 1
}

########################################################################################################################
PREFIX=""
while [ $# -ge 1 ]; do
    case $1 in
        -h|-\?|--help)
            show_help
            exit 0
            ;;
        -p|--prefix)
            if [ $# -ge 2 ]; then
                PREFIX="$2"
                shift 2
            else
                throw_error "--prefix requires a non-empty argument"
            fi
            ;;
        --PREFIX=?*)
            PREFIX=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --)   # explicit call to end of all options
            shift
            break
            ;;
        -?*)  # unsupported option
            throw_error "Unknown option \"$1\". use --help for syntax"
            ;;
        *)  # not an option, a positional argument. break out
            break
            ;;
    esac
done

if [[ "$#" -lt 5 ]]; then
    show_help && exit 0
fi

########################################################################################################################
mer_sz=$1

hap_a_db=$2    # meryl hap_a database dir
hap_b_db=$3    # meryl hap_b database dir

a_fa=$4     # hap_a.fasta/fa assembly
b_fa=$5     # hap_b.fasta/fa assembly
########################################################################################################################
a_name=$(basename "${a_fa}" | sed 's/.fasta$//' | sed 's/.fa$//')
b_name=$(basename "${b_fa}" | sed 's/.fasta$//' | sed 's/.fa$//')
tmp=$(basename -- "${hap_a_db}")
hap_a_name="${tmp%.*}"
tmp=$(basename -- "${hap_b_db}")
hap_b_name="${tmp%.*}"

num_cores=$(grep -c ^processor /proc/cpuinfo)
machine_mem=$(echo $(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024))))
meryl_mem=$(($machine_mem - 2))

echo "=============================="
echo "Start counting"
if [[ ! -f "${PREFIX}_${a_name}_${hap_a_name}.counts" ]]; then
    echo "${a_name} in ${hap_a_name}"
    date -u
    meryl-lookup -existence -threads "${num_cores}" -memory "${meryl_mem}" \
        -min "${mer_sz}" \
        -max "${mer_sz}" \
        -mers "${hap_a_db}" \
        -sequence "${a_fa}" | \
        awk -v name="${a_name}" 'BEGIN{OFS="\t"} {print name, $0}' \
        > "${PREFIX}_${a_name}_${hap_a_name}.counts"
    date -u
    echo
fi

if [[ ! -f "${PREFIX}_${b_name}_${hap_a_name}.counts" ]]; then
    echo "${b_name} in ${hap_a_name}"
    date -u
    meryl-lookup -existence -threads "${num_cores}" -memory "${meryl_mem}" \
        -min "${mer_sz}" \
        -max "${mer_sz}" \
        -mers "${hap_a_db}" \
        -sequence "${b_fa}" | \
        awk -v name="${b_name}" 'BEGIN{OFS="\t"} {print name, $0}' \
        > "${PREFIX}_${b_name}_${hap_a_name}.counts"
    date -u
    echo
fi

if [[ ! -f "${PREFIX}_${a_name}_${hap_b_name}.counts" ]]; then
    echo "${a_name} in ${hap_b_name}"
    date -u
    meryl-lookup -existence -threads "${num_cores}" -memory "${meryl_mem}" \
        -min "${mer_sz}" \
        -max "${mer_sz}" \
        -mers "${hap_b_db}" \
        -sequence "${a_fa}" | \
        awk -v name="${a_name}" 'BEGIN{OFS="\t"} {print name, $0}' \
        > "${PREFIX}_${a_name}_${hap_b_name}.counts"
    date -u
    echo
fi

if [[ ! -f "${PREFIX}_${b_name}_${hap_b_name}.counts" ]]; then
    echo "${b_name} in ${hap_b_name}"
    date -u
    meryl-lookup -existence -threads "${num_cores}" -memory "${meryl_mem}" \
        -min "${mer_sz}" \
        -max "${mer_sz}" \
        -mers "${hap_b_db}" \
        -sequence "${b_fa}" | \
        awk -v name="${b_name}" 'BEGIN{OFS="\t"} {print name, $0}' \
        > "${PREFIX}_${b_name}_${hap_b_name}.counts"
    date -u
    echo
fi

echo "=============================="
echo "Merge"

paste "${PREFIX}_${a_name}_${hap_a_name}.counts" "${PREFIX}_${a_name}_${hap_b_name}.counts" > "${PREFIX}_${a_name}.counts"
paste "${PREFIX}_${b_name}_${hap_a_name}.counts" "${PREFIX}_${b_name}_${hap_b_name}.counts" > "${PREFIX}_${b_name}.counts"

echo -e "Assembly\tContig\t${hap_a_name}\t{$hap_b_name}\tTotal" > "${PREFIX}_hapmers.counts"
cat "${PREFIX}_${a_name}.counts" "${PREFIX}_${b_name}.counts" | \
    awk 'BEGIN{OFS="\t"} {print $1, $2, $5, $NF, $3}' >> "${PREFIX}_hapmers.counts"

echo "=============================="
echo "             DONE             "
echo "=============================="
