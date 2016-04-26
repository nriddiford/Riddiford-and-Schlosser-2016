################
### de.genes ###
################

#!/usr/bin/perl
use strict; 
use warnings;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

# This script will attach differenital expression info from the cufdiff output file 'gene_exp.diff' to transcripts assembled in 'transcripts.pl'

# Run as 'perl de.genes.pl' 


open my $in, '<', 'diff_out/gene_exp.diff' or die "usage perl de.genes.pl\nIs the directory 'diff_out'in the current directory?\n";

my %diff_out;
while (<$in>){
	chomp;
	next if /test_id/;
	my @split = split(/\t/);
	my $xloc = $split[0];
	my $scaffold = $split[3];
	my ($control_val, $expt_val, $change_val)  = @split[7,8,9];
	$diff_out{$xloc} = [$scaffold, $control_val, $expt_val, $change_val];
}

open my $hits, '<', 'transcripts.txt'  or die "'transcripts.txt' must be in the current directory\n";

my (%genes, $xloc);
while (<$hits>) {
	chomp;
	($xloc) = /(XLOC_\d+)/ if /^>/;  
    $genes{$xloc} = $_ if /^[A-Z]/;	
}

open my $de_genes_print, '>', 'de_transcripts_new.txt' or die "$!\n";

foreach my $key (sort keys %diff_out){
	my ($scaffold, $control_val, $expt_val, $change_val) = @{$diff_out{$key}};
	print $de_genes_print ">$key\_$control_val:$expt_val,_Change:$change_val\n$genes{$key}\n" if exists $genes{$xloc};
}