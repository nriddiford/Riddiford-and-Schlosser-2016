#################
### Pearson's ###
#################
 
#!/usr/bin/perl
use warnings;
use strict; 
use feature qw(say);
use Data::Dump 'dump';
use Data::Dumper;

# This script will read in fpkm values from cuffdiff output file 'genes.read_group_tracking' and open a .csv file containing:
# FPKM rep 1	FPKM rep 2
# This can then be used as input for a standard Pearson's correlation (I did this in R) 

# Run as perl pearsons.pl 
# Will open output in excel

open my $in, '<', "/Users/Nick/Documents/Bioinformatics/Data/eya1/diff_out/genes.read_group_tracking" or die "usage perl pearsons.pl\n";

my %data;
while (<$in>){
	chomp;
	next unless /^XLOC/;
	my @split = split(/\t/);
	$data{$split[0]}{$split[1]}{$split[2]} = $split[6];
}

open my $out, '>', 'output.csv';

print $out "fpkm1,fpkm2\n";

my (%rep_count, $log_trans_fpkm, %filter, $fpkm);

for my $xloc (sort keys %data){
	for my $condition (keys $data{$xloc}){
		for my $replicate (sort keys $data{$xloc}{$condition}){
			$fpkm = $data{$xloc}{$condition}{$replicate};
			
			next if $fpkm == 0;
			$rep_count{$xloc}{$condition}++;
			
			$filter{$xloc}{$condition}{$replicate} = $fpkm if $replicate == 0 and $rep_count{$xloc}{$condition} == 1 and $condition eq 'Echx';
			$filter{$xloc}{$condition}{$replicate} = $fpkm if $replicate == 1 and $rep_count{$xloc}{$condition} == 2 and $condition eq 'Echx';
		}
		delete $filter{$xloc} if $fpkm == 0 or $fpkm > 5000;
	}
}

for my $id (sort keys %filter){
	for my $condition (keys $filter{$id}){
		for my $replicate (sort keys $filter{$id}{$condition}){
			print $out "$filter{$id}{$condition}{$replicate}," if $replicate == 0 and $condition eq 'Echx';
			print $out "$filter{$id}{$condition}{$replicate}" if $replicate == 1 and $condition eq 'Echx';
		}
		print $out "\n";
	}		
}

system("open output.csv");
