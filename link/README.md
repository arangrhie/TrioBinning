# Expanded trio binning though k-mer linking

This README contains a step-by-step instruction to do the expanded trio-binning on HiFi reads.
Note that this is a pseudocode for the entire procedure, and has only being tested on HG002 trio data from the T2T/GIAB.

## Requires
* Meryl v1.4
* Merqury, latest tip
* Canu for `splitHaplotype` (tested on r10354 790082e884451f5af19831ed50eb661143a9a300)
* GenomeScope2

## 1. Prepare k-mer dbs
### 1. Collect 21-mers of the parents
```shell
meryl count k=21 mat.*.fastq.gz output mat.meryl
meryl count k=21 pat.*.fastq.gz output pat.meryl
meryl count k=21 child.*.fastq.gz output child.meryl
```

### 2. Collect hapmers
```shell
$MERQURY/trio/hapmers.sh mat.meryl pat.meryl child.meryl
```
Output files needed:
* mat.hapmer.meryl
* pat.hapmer.meryl
* inherited_hapmers.png

### 3. GenomeScope2
```shell
meryl histogram child.meryl > child.hist
Rscript $tools/genomescope2.0/genomescope.R \
  -i child.hist -o sample_gs2 -p 2 -k 21 --fitted_hist
```
Output files needed:
* lookup_table.txt
* fitted_hist.png

### 4. Collect 1-copy region
```shell
LO=`awk -F "," '$1==1 {print (NR-1); exit}' lookup_table.txt`
HI=`awk -F "," '$1==2 {print NR; exit}'     lookup_table.txt`
echo "LO: $LO, HI: $HI"
```
Multiplicity of the 1-copy region for the child X will be
* `$LO` < X < `$HI`

This is the multiplicity range that GenomeScope2 predicts the chance for a k-mer belonging to the 1-copy region
is higher than belonging to an erroneous or 2-copy region of the genome.

## 2. `splitHaplotype` using hapmers

```shell
/path/to/canu/build/bin/splitHaplotype \
 -cl 1000 -memory 48 -threads 24 \
 -H ../mat.hapmer.meryl 1 mat.fasta.gz \
 -H ../pat.hapmer.meryl 1 pat.fasta.gz \
 -A ./unk.fasta.gz \
 `cat input.fofn | awk '{print "-R "$0}' | tr '\n' ' '`
 ```
 * input.fofn: HiFi reads of the child
 * `-cl 1000`: Reads <1kb will not be reported.
 
## 3. Link hapmers
### 1. Collect k-mers from the binned reads
This time, we are collecting k-mers of the binned HiFi reads.
```shell
meryl count k=21 mat.fasta.gz   output read-mat.meryl
meryl count k=21 pat.fasta.gz   output read-pat.meryl
```

### 2. Link hapmers from the reads
This step will be repeated at least twice.
Let's do this in a clean dir so we don't mix up with results from previous rounds.

```shell
mkdir round1
cd round1

ln -s /path/to/mat.hapmers.meryl
ln -s /path/to/pat.hapmers.meryl
ln -s /path/to/child.meryl
```

Check `HI` and `LO` are still valid from [step 1-4](https://github.com/arangrhie/TrioBinning/new/master#4-collect-1-copy-region).
```shell
echo -e "mat-lnk $LO" >  cutoffs.txt
echo -e "mat-lnk $HI" >> cutoffs.txt
cat cutoffs.txt
```

Let's collect the linked hapmers.
```shell
link_hapmers() {
  hap=$1
  other_hap=$2

  echo "Processing $hap"

  # hapmers in 1-copy range from the initial binned reads
  meryl output read-$hap.1cp.meryl \
  [ less-than $HI [ greater-than $LO read-$hap.meryl ] ]

  # exclude hapmers from both haplotypes to get the newly added hapmers
  meryl difference \
    [ difference read-$hap.1cp.meryl $hap.hapmer.meryl ]\
    $other_hap.hapmer.meryl \
    output read-$hap.1cp.only.meryl

  # union the hapmers and newly binned mers
  meryl union-sum $hap.hapmer.meryl read-$hap.1cp.only.meryl \
    output $hap.lnk.meryl
}

link_hapmers mat pat
link_hapmers pat mat

meryl union-sum  mat.lnk.meryl pat.lnk.meryl     output mat_pat.lnk.meryl
meryl difference child.meryl   mat_pat.lnk.meryl output child_only.meryl
```

### 3. Plot linked hapmers
Plot to double check the results look reasonable.
```shell
# Obtain .hist for spectra-cn.R
$MERQURY/plot/to_hist_for_plotting.sh \
  child_only.meryl unbinned \
  mat.lnk.meryl mat-lnk \
  pat.lnk.meryl pat-lnk \
> linked_hapmers.hist

Rscript $MERQURY/plot/plot_spectra_cn.R -f linked_hapmers.hist \
 -o linked_hapmers \
 -m 150 \
 -n 65000000 \
 -l cutoffs.txt
 ```
Output:
* linked_hapmers.png

Note that `-m` and `-n` are adjusted for plotting. Try first without it, and choose the appropriate X and Y axis maximum to plot.

Compare `linked_hapmers.png` with `inherited_hapmers.png` from [step 1-2](https://github.com/arangrhie/TrioBinning/new/master#2-collect-hapmers).

### 4. Re-bin the `unk.fasta.gz` from previous the round
```shell
# Collect input.fofn to run splitHaplotype
ls /path/to/unk.fasta.gz > input.fofn
```
Now run `splitHaplotype` as in [step 2](https://github.com/arangrhie/TrioBinning/new/master#2-splithaplotype-using-hapmers).

## 4. Linking, round 2
### 1. Collect k-mers from the newly binned reads
```shell
mkdir round2
cd round2

ln -s /path/to/child.meryl

for hap in mat pat
do  
  # Collect k-mers from the newly binned reads
  meryl count k=21 ../$hap.fasta.gz output read-$hap.meryl

  # Add them to the previous ones
  meryl union-sum read-$hap.meryl ../../read-$hap.meryl \
    output read-$hap.all.meryl
done
```
Output needed:
* read-mat.all.meryl
* read-pat.all.meryl

### 2. Expand linked k-mers
```shell
linkend_hapmers() {
  hap=$1
  other_hap=$2

  READ_HAP=read-$hap.all.hapmer.meryl
  READ_HAP_1CP=read-$hap.1cp.meryl
  READ_HAP_1CP_NEW=read-$hap.1cp.only.meryl
  HAP_LNK=$hap.lnk.meryl
  HAP_LNK_OTHER=$other_hap.lnk.meryl

  echo "Processing $hap"

  # hapmers in 1-copy range from the initial binned reads
  meryl output $READ_HAP_1CP \
  [ less-than $HI [ greater-than $LO $READ_HAP ] ]

  # exclude linked hapmers from the other haplotype of the previous round to get the newly added hapmers
  meryl difference \
    [ difference $READ_HAP_1CP ../$HAP_LNK ]\
    ../$HAP_LNK_OTHER \
    output $READ_HAP_1CP_NEW

  # union the linked hapmers from prev. round and newly binned mers
  meryl union-sum ../$HAP_LNK $READ_HAP_1CP_NEW \
    output $HAP_LNK
}

linkend_hapmers mat pat
linkend_hapmers pat mat

meryl union-sum mat.lnk.meryl pat.lnk.meryl output mat_pat.lnk.meryl
meryl difference child.meryl mat_pat.lnk.meryl output child_only.meryl
```
Output needed:
* mat.lnk.meryl
* pat.lnk.meryl
* child_only.meryl

### 3. Plot and re-bin
Use the `cutoff.txt` file from previous round, and run the same plotting script and binning as in [this](https://github.com/arangrhie/TrioBinning/new/master#3-plot-linked-hapmers)
and [this](https://github.com/arangrhie/TrioBinning/new/master#4-re-bin-the-unkfastagz-from-previous-the-round) steps.

## 5. Evaluate binning and repeat
Pick up `LO` and `HI` from [here](https://github.com/arangrhie/TrioBinning/new/master#4-collect-1-copy-region) or `cutoffs.txt`.
`PEAK` is the multiplicity of the maximum counts in the 1-copy kmers.
Obtain histogram for `$hap.lnk.meryl` with `meryl histogram`, and pick the value in 1st column (multiplicity, between `LO` and `HI`) where the 2nd column peaks. 

```shell
OUT=`meryl statistics  child.meryl | head -n3 | tail -n1 | awk '{print $2}'`
echo "child.meryl ALL: $OUT"

OUT=`meryl histogram  child.meryl \
  | awk -F "\t" -v LO=$LO -v HI=$HI '{if (LO<$1 && $1<HI) sum+=$2} END {print sum}'`
echo "child.meryl in 1-copy region ($LO, $HI): $OUT"

OUT=`meryl histogram  child.meryl \
  | awk -F "\t" -v PK=$PEAK '$1==PK {print $2}'`
echo "child.meryl at $PEAK: $OUT"

# Repeat this on the initial hapmers as well
for hap in mat pat
do
  OUT=`meryl statistics  $hap.lnk.meryl | head -n3 | tail -n1 | awk '{print $2}'`
  echo "$hap.lnk.meryl ALL: $OUT"

  OUT=`meryl histogram  $hap.lnk.meryl \
  | awk -F "\t" -v LO=$LO -v HI=$HI '{if (LO<$1 && $1<HI) sum+=$2} END {print sum}'`
  echo "$hap.lnk.meryl in 1-copy region ($LO, $HI): $OUT"

  OUT=`meryl histogram  $hap.lnk.meryl \
  | awk -F "\t" -v PK=$PEAK '$1==PK {print $2}'`
  echo "$hap.lnk.meryl at $PEAK: $OUT"
done
```

Compare values of linked hapmers (`$hap.lnk.meryl`) at each round against the values from `child.meryl`.
Propotion for ( mat.lnk.meryl + pat.lnk.meryl / child.meryl ) should increase at each round compared to the values from the initial hapmers.

_NOTE_: Keep in mind that the 1-copy region contains k-mers from the erroneous and 2-copy region, the proportion will never reach 100%.

The end of the log file from `splitHaplotype` shows stats of the binned reads (num. reads and num. bases).
Decide if you want to repeat [step 4](https://github.com/arangrhie/TrioBinning/new/master#4-linking-round-2) given the increase.
