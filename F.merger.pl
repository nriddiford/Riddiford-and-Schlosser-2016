##############
### Merger ###
##############

#!/usr/bin/perl
use warnings;
use strict; 
# use diagnostics;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use feature qw(say);
use autodie;

# This script will annotate 'gene_exp.diff' information on the direction of change
# To be used as input for scatterplot

# Run as perl merger.pl 

open my $in, '<', 'gene_exp.diff';

my (%data, $sig);
while (<$in>){
	chomp;
	next unless /^XLOC/;
	my @split = split(/\t/);
	$data{$split[0]} = $split[13];	
}

open my $fpkm, '<', 'new.txt';

open my $out, '>', 'merged.txt';

open my $genes, '<', '/Users/Nick/Documents/Bioinformatics/Data/Merged_reps/Blast/merge.merged.blast_out.20.1.15.txt';

while (<$genes>) { 
	chomp;
	next if /^\s*$/;
	my ($xloc) = /(XLOC_\d+)/;
	my ($change) = /Change:(-?\d+\.\d+|-?inf|-?\d+)/;
	my ($q) = /q:(\d+\.\d+|-?\d+)\s/;
	my ($container) = /XLOC_\d+_(.+:.+),_Change/g;
	my ($con_val, $expt_val) = split(/:/, $container);
	my @split = split(/\t/);
	my ($condition, $percent_id, $gene) = ($split[0], $split[2], $split[3]);
	}	
}

print $gene_names "gene\tcontrol\texpt\n";
while (<$fpkm>){
	chomp;
	my @split = split(/\t/);
	my ($con, $exp) = ($split[9], $split[13]);
	if ($data{$split[0]}) {
		if ($data{$split[0]} eq 'yes' and $con < $exp){
			print $out "$_\t$data{$split[0]}\tup\n";
		}
		elsif ($data{$split[0]} eq 'yes' and $con > $exp){
 			print $out "$_\t$data{$split[0]}\tdown\n";
		}
		elsif ($data{$split[0]} eq 'no'){
			print $out "$_\t$data{$split[0]}\tno\n";
		}
	}
}
