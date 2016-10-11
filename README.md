# Bioinformatic pipeline for Riddiford and Schlosser 2016

Please cite [Riddiford and Schlosser 2016](https://elifesciences.org/content/5/e17666) when you use this pipeline for data analysis

# Table of Contents

* [Trimming and mapping](#trimming-and-mapping)
  * [FastQC](#run-fastqc-to-visually-inspect-all-sequencing-results)
  * [Trimmomatic](#run-trimmomatic-to-quality-filter-reads)
  * [Tophat](#run-tophat2-on-each-set-of-paired-reads)
* [Transcript assembly and differential expression analysis](#transcript-assembly-and-differential-expression-analysis)
  * [Cufflinks](#assemble-transcripts-and-calcuulate-abundance-estimation-with-cufflinks)
  * [Cuffmerge](#merge-assemblies-using-cuffmerge)
  * [Cuffdiff](#differential-expression-using-cuffdiff)
* [Estimating variance between biological replicates](#estimating-variance-between-biological-replicates)
* [Annotating transcript models](#annotating-transcript-models)
  * [Build transcript models](#build-transcript-models)
  * [Add differential expression information](#add-differential-expression-information)
  * [Annotate](#annotate)
* [Gene Enrichment Analysis](#gene-enrichment-analysis)
* [Find co-differentially expressed genes](#find-co-differentially-expressed-genes)
  * [For single conditions](#for-single-conditions)
  * [For merged conditions](#for-merged-conditions)
* [Gene Ontology](#gene-ontology-analysis-on-discrete-gene-sets)
  * [DAVID](#convert-to-human-orthologs)


## Trimming and mapping 


### Run fastqc to visually inspect all sequencing results

Run from the same directory as fasta files. Will output .fastq files with same name

```{perl}
#!/usr/bin/perl
use warnings;
use strict;
use feature qw(say);

open (FILES, "ls *.txt |");
while (<FILES>) {
	chomp;
	my ($file) = (split)[0];
	say "Parsing $file...";
	system ("fastq_quality_trimmer -t 30 -l 75 -i $file -o $file.fastq");
}
```


### Run Trimmomatic to quality filter reads

`filter_set` contains list of primer sequences to exclude

```{java}
java -classpath /path/to/Trimmomatic/trimmomatic-0.25.jar org.usadellab.trimmomatic.TrimmomaticPE
-threads 12 \
-phred33 \
<pe_1> <pe_2> <paired_output_1> <unpaired_output_1> <paired_output_2> <unpaired_output_2> \
ILLUMINACLIP:<filter_set> \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:12 MINLEN:36
```

Build bowtie index

Run in same directory as `genome.fasta` file
```
bowtie2-build <genome.fasta>
```


### Run Tophat2 on each set of paired reads

``` tophat2 -p 12 -r 250 -N 3 -o <ouput_dir> <bowtie-index> <pe_1> <pe_2>```

Get stats for run and rename .bam file

```cd <output_dir>```

```samtools flagstat accepted_hits.bam > stats.txt```

```mv accepted_hits.bam <condition.rep.bam>```


## Transcript assembly and differential expression analysis


### Assemble transcripts and calcuulate abundance estimation with Cufflinks

```cufflinks -p 12 -b <bowtie-index> -o <ouput_dir> <condition.rep.bam>```


### Merge assemblies using Cuffmerge

```ls -1 <path/to/cufflinks_output_condition_1/transcripts.gtf> <condition_n/transcripts.gtf> > assemblies.txt ```
```cuffmerge assemblies.txt```


### Differential expression using Cuffdiff

```
cd <cufflinks_output_dir> \
cuffdiff -p 12 \
-o <cufflinks_output_dir/diff_out> \
-b <bowtie-index/.fasta> \
-L <L1,L2> <merged.gtf> <L1_rep1.bam>,<L1_rep2.bam> <L2_rep1.bam>,<L2_rep2.bam>
````


## Estimating variance between biological replicates

Run **pearsons.pl** from same directory as Cuffdiff output file `genes.read_group_tracking`

```perl pearsons.pl```

Output can used as input for a standard Pearson's correlation (e.g. in R)


## Annotating transcript models


### Build transcript models

To build transcript models from genome using cufflinks output, run **transcripts.pl** from same directory as `merged.gtf`

```perl transcripts.pl```

This script will read in exonic positions from cufflinks output file `merged.gtf` and relate then to genomic positions. The program outputs:

 * All transcript variants for a given gene
 * All exons in each gene
 * The longest transcript variant for all assembled genes (numbered by transcript variant)


### Add differential expression information

Run **de.genes.pl**. Will need to moify file locations in script

```perl de.genes.pl```

This script will attach differenital expression info from the cufdiff output file `gene_exp.diff` to transcripts assembled in **transcripts.pl**


### Annotate

BLAST output from **de_genes.pl** (`de_transcripts_.fa`) against Xenopus mRNA database:

```
blastn -db <path_to_Xenopus_DB> \
-query <de_transcripts_.fa> \
-num_threads 12 \ 
-evalue '0.00001' \
-perc_identity 80 \
-outfmt "6 qseqid pident sseqid" \
-task blastn \
-max_target_seqs 1 | sort -u -k1,1 > <output.txt>
```

Run **annotator.pl**

```perl annotator.pl <Xenopus blast output> <condition> <de_genes.txt>```


## Gene Enrichment Analysis

To look for gene enrichement between two conditions, run **enrichment.pl**

```perl enrichment.pl <blast.file1> <blast.file2>```

This script will read in differentially expressed genes (output from **blast.pl**) for two conditions and output the top 10% for each list, as well as input data for a Chi-squared table 

Perform chi squared test using output from **enrichment.pl** (e.g. use http://www.socscistatistics.com/tests/chisquare/)


## Find co-differentially expressed genes 


### For single conditions

Run **godzilla.pl** for blast output for multiple conditions. 
"Control" must always be listed first

Blast output must have the following, tab delimited output:
`gene	XLOC_control.fpkm:expt.fpkm,_Change:FC_val	%id	blast hit`

```perl godzilla.pl control.txt six.txt six-eya.txt eya.txt```

This program will take N blast hit results and output those genes which are found to be overlapping in `$n` conditions.


### For merged conditions

Run **DE_sig.pl** 

"Control" must always be listed first

Blast output must have the following, tab delimited output:
`gene	XLOC_control.fpkm:expt.fpkm,_Change:FC_val	%id	blast hit`

```perl DE_sig.pl control.txt six.txt six-eya.txt eya.txt```


# Gene ontology analysis on discrete gene sets

Run **sets.pl** to look for common genes between different de gene sets (lists outputted from **godzilla.pl**)

```perl sets.pl list1 list2 list3```


# Convert to Human orthologs

Blast output (e.g. 'setA') against Human Uniprot DB

```
blastx -db <path_to_Human_UniDB> \
-query <setA.fa> \
-num_threads 12 \
-evalue '0.001' \
-outfmt "6 qseqid pident sseqid" \
-max_target_seqs 1 | sort -u -k1,1 > <output.txt>
```

Resulting output can then be used as input for [David](https://david.ncifcrf.gov)
