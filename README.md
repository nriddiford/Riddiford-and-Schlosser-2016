# Pipeline for Riddiford and Schlosser 2016

Please cite [Riddiford and Schlosser 2016](https://elifesciences.org/content/5/e17666) when you use this pipeline for data analysis

# Trimming and mapping 

Run **fastqc** to visually inspect all sequencing results

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

Run **Trimmomatic** to quality filter reads

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

Run **Tophat2** on each set of paired reads


``` tophat2 -p 12 -r 250 -N 3 -o <ouput_dir> <bowtie-index> <pe_1> <pe_2>```

Get stats for run and rename .bam file

```cd <output_dir>```

```samtools flagstat accepted_hits.bam > stats.txt```

```mv accepted_hits.bam <condition.rep.bam>```

# Transcript assembly and differential expression analysis

Assemble transcripts and calcuulate abundance estimation with **Cufflinks**

```cufflinks -p 12 -b <bowtie-index> -o <ouput_dir> <condition.rep.bam>```

Merge assemblies using **Cuffmerge**

```ls -1 <path/to/cufflinks_output_condition_1/transcripts.gtf> <condition_n/transcripts.gtf> > assemblies.txt ```
```cuffmerge assemblies.txt```

Differential expression using **Cuffdiff**

```
cd <cufflinks_output_dir> \
cuffdiff -p 12 \
-o <cufflinks_output_dir/diff_out> \
-b <bowtie-index/.fasta> \
-L <L1,L2> <merged.gtf> <L1_rep1.bam>,<L1_rep2.bam> <L2_rep1.bam>,<L2_rep2.bam>
````

# Estimating variance between biological replicates

Run **pearsons.pl** from same directory as Cuffdiff output file `genes.read_group_tracking`

Output can used as input for a standard Pearson's correlation (e.g. in R)

# Annotating transcript models

To build transcript models from genome using cufflinks output, run **transcripts.pl** from same directory as `merged.gtf`

This script will read in exonic positions from cufflinks output file `merged.gtf` and relate then to genomic positions. The program outputs:

 * All transcript variants for a given gene
 * All exons in each gene
 * The longest transcript variant for all assembled genes (numbered by transcript variant)


Run **de.genes.pl**. Will need to moify file locations in script

This script will attach differenital expression info from the cufdiff output file `gene_exp.diff` to transcripts assembled in **transcripts.pl**


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

Run annotator.pl

# 10. Gene Enrichment Analysis 
# run enrichment.pl
# perform chi squared test using output from 'enrichment.pl' ( e.g. use http://www.socscistatistics.com/tests/chisquare/)

# 11. Finding co-differentially expressed genes (for single conditions)
# run godzilla.pl

11. Finding co-differentially expressed genes (for merged conditions)
# run DE_sig.pl

# 12. Gene ontology analysis on discrete gene sets
# run sets.pl
# Blast output (e.g. 'setA') against Human Uniprot DB

blastx -db <path_to_Human_UniDB> \
-query <setA.fa> \
-num_threads 12 \
-evalue '0.001' \
-outfmt "6 qseqid pident sseqid" \
-max_target_seqs 1 | sort -u -k1,1 > <output.txt>
