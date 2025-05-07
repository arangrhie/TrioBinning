# TrioBinning
Scripts for manually classifying reads and make figures

## Dependencies
* meryl and meryl-lookup from [Meryl v1.4.1](https://github.com/marbl/meryl)
* hapmer dbs built with Merqury [hapmers.sh](https://github.com/marbl/merqury/blob/master/trio/hapmers.sh)
* Set env variable for tools and meryl in `~/.bash_profile`
	```sh
	export tools=/path/to/tools # where all repos are
    export PATH=$PATH:$tools/meryl/bin
	```

## Assign
* Assign each read to a haplotype using  `assign.sh`.
The script runs `meryl-lookup -existence -sequence $fq -mers $hapmer1 $hapmer2` underneath and prints a tsv file with the following:
* normalized count for hap1 (`$4/$3`), normalized count for hap2 (`$6/$5`)
* set `score1` to be always larger or equal to `score2`, `hap` to the corresponding haplotype matching the scores
* if `score2` is 0, `score1` is larger than 0, assign to `hap1`
* if `score2` is larger than 0, and `score1/score2` > `cutoff`, assign to `hap2`
* otherwise, assign to `unknown`

Submit via slurm array job to run on each line:
```sh
cutoff=4.0 # optional, default is 4.0
./_submit_assign.sh mat.hapmer.meryl pat.hapmer.meryl input.fofn $cutoff
```

## Plot hapmers in bubbles
Visually inspect how the reads are assigned by using the bubble plot, adjust `$cutoff` if needed.
```sh
cat $input.assigned.*.$cutoff > reads_assigned.$cutoff
Rscript hapmer_plot_alt.R --file reads_assigned.$cutoff --output reads_assigned.$cutoff.png
```

## Bin reads
Bin input fq reads according to the assigned.$cutoff file. Submit via slurm array job to run on each line.
```sh
./_submit_binning.sh input.fofn $cutoff
```
