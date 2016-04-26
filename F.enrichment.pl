##################
### Enrichment ###
##################

#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use feature 'say';
use POSIX;

# This script will read in differentially expressed genes (output from 'blast.pl') for two conditions
# and output the top 10% for each list, as well as input data for a Chi-squared table 

# Run as perl enrichment.pl <blast.file1> <blast.file2>


if ($#ARGV < 0) {
	die "usage enrichment.pl <blast.file1> <blast.file2>";
	exit;
}

my ($g1, $g2);
my $p = 10;
my ($r1, $r2);
my (%top1, %top2);

my $threshold = 0.5;

my (%duplicates, %gene_list, %change);
my 	$fileno;

my $all_genes;

blast_extractor($_) foreach @ARGV;

sub blast_extractor {
    my $file = shift;
    open my $output, '<', $file or die "Can't read file '$file' [$!]\n";
	$fileno++;
	my @change;
	while (<$output>) {
		chomp;
		next if /^\s*$/;
		my ($xloc) = /(XLOC_\d+)/;
		my ($change) = /Change:(-?\d+\.\d+|-?inf|-?\d+)/;
		my ($container) = /XLOC_\d+_(.+:.+),_Change/g;
		my ($con_val, $expt_val) = split(/:/, $container);
		my @split = split(/\t/);
		my ($condition, $percent_id, $gene) = ($split[0], $split[2], $split[3]);
		$gene =~ s/_/ /g;
		$gene =~ s/,/ /g;
		my @small_split = split(/\|/, $gene);
		my $small_gene = $small_split[4];
		my ($count) = $duplicates{$small_gene}++;
		$gene_list{$small_gene}{$count} = [$xloc, $con_val, $expt_val, $change];
	}

	my (%filtered, %expt_count, %pos_count, %neg_count, %con_sum, %expt_sum, %genes_by_count, %all, @xlocs);
	my $change_sum;

	my $gene;
	
	foreach $gene (sort keys %gene_list) {
			foreach my $count (sort keys $gene_list{$gene} ){
				my ($xloc, $con_val, $expt_val, $change_val) = @{$gene_list{$gene}{$count}};
				# Remove the occurrences of 'inf' log(2) values by adding 0.001 to each raw value of 0
				$con_val = ($con_val + 0.001) if $con_val == 0;
				$expt_val = ($expt_val + 0.001) if $expt_val == 0;

				# Foreach gene/condition sum up all the raw values
				$con_sum{$gene} += $con_val;
				$expt_sum{$gene} += $expt_val;

				$con_val = $con_sum{$gene};
	 			$expt_val = $expt_sum{$gene};

				my $raw_val = $expt_val/$con_val;
		 		$change_sum = eval sprintf('%.6f', (log($raw_val)/ log(2)));

				$genes_by_count{$gene} = [$xloc, $con_val, $expt_val, $change_sum, $count];

	 			delete $genes_by_count{$gene} if abs $change_sum < $threshold or $genes_by_count{$gene}[2] < 1;
			}
	}

	for my $g (sort keys %genes_by_count) {	
		my ($x, $c, $e, $fc, $co) = @{$genes_by_count{$g}};
		push @change, [$g, $fc];
	}
	
	my @top;
	my $top10;
	my $counter = 0;
	
	for (sort { $b->[1] <=> $a->[1] } @change){
		$counter++;
		$top10 = int(10/100 * ($#change));
		push @top, @$_[0] if $counter < $top10;
	}	
	
	
	$r1 = $top10 if $fileno == 1;
	$r2 = $top10 if $fileno == 2;
	$g1 = $#top if $fileno == 1;
	$g2 = $#top if $fileno == 2;
	
	say "$top10 genes in top $p% of $file";
	say "$counter genes in $file";	
	
	$all_genes += $counter;
	$all_genes -= $top10;
	$top10 = ();
	

	if ($fileno == 1) {$top1{$_} = 1 foreach @top}
	if ($fileno == 2) {$top2{$_} = 2 foreach @top}

@change=();

delete $genes_by_count{$_} for keys %genes_by_count;
delete $gene_list{$_} for keys %gene_list;
}

my $overlap = 0;
my (@overlap, %overlap);
foreach my $gene (keys %top1){
	if (exists $top2{$gene}){
		$overlap++;
		push @overlap, $gene;
		$overlap{$gene} = 1;
	}  		
}

say "\n'a' -> Genes in top $p% for condition1 & condition2 = $overlap";
my $not_overlap1 = ceil($r1 - $overlap);
my $not_overlap2 = ceil($r2 - $overlap);

my $non_overlap = ceil(($g1 + $g2)/2);
say "'b' -> Genes in the top $p% for condition2, but not in the top $p% for condition1  = $not_overlap1";
say "'c' -> Genes in the top $p% for condition1, but not in the top $p% for condition2  = $not_overlap2";
say "'d' -> Genes in neither the top $p% for condion1 nor condition2 = $all_genes ";


print "-----------------------\n";
print "|          |           |\n";
print "|    $overlap     |    $not_overlap1    |\n";
print "|          |           |\n";
print "-----------------------\n";
print "|          |           |\n";
print "|    $not_overlap2   |   $all_genes    |\n";
print "|          |           |\n";
print "-----------------------\n";


print "\n";
