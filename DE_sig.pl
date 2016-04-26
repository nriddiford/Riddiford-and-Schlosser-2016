#!/usr/bin/perl
use warnings;
use strict; 
# use diagnostics;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use feature qw(say);
use autodie;

##############
#### INTRO ###
##############

# This program will take N blast hit results and output those genes which are found to be overlapping in $n conditions.
# The test data be run as follows: perl overlapping.pl control.txt six.txt six-eya.txt eya.txt

# This progeram can be used to generate various gene lists:

### Downreg ###
# line 108 'next if exists $neg_count{$gene};' => 'next if exists $pos_count{$gene};'
# line 249 '$change_val = 0 if $change_val < 0;' => '$change_val = 0 if $change_val > 0;'

### Single gene lists ###
# Only input control vs condition of interest
# line 48 => $n = 1
# line 243 'next unless $seen{$gene} >= $n;' => 'next unless $seen{$gene} == $n;'

### Intersection between only 2 condtitions ### 
# line 48 => $n = 2
# line 243 'next unless $seen{$gene} >= $n;' => 'next unless $seen{$gene} == $n;'
# line 59 'my $filter = 'six1-eya1';' => enable this option, and specify the non-intersecting condition
# line 111 'next if $genes_by_count{$gene}{$filter};' => enable this option

##########
##########

if ($#ARGV < 0) {
	usage();
	exit;
}

my (%duplicates, %gene_list, %change, %seen, %seen_all);

print "Opening blast ouput foreach condition...\n";

blast_extractor($_) foreach @ARGV;

print "Enter direction: ";
chomp(my $direction = <STDIN>);

my $n = 1;

# print "Enter query: ";
# chomp(my $query = <STDIN>);
# die "no search string entered\n" unless $query =~ /^[A-z]/;


my (%filtered, %expt_count, %pos_count, %neg_count, %con_sum, %expt_sum, %genes_by_count, %all, @xlocs);

my $change_sum;

# my $filter = 'six1-eya1';
my $threshold = 1;

say "Filtering genes for DE genes...";

# say "Returning matches for $query";

my %qs;
for my $gene ( keys %gene_list) {	
	foreach my $con (keys $gene_list{$gene}) {
		$expt_count{$gene}++ unless $con eq 'control';
		foreach my $count (sort keys $gene_list{$gene}{$con} ){
			my ($xloc, $con_val, $expt_val, $change_val, $q) = @{$gene_list{$gene}{$con}{$count}};
			# Remove the occurrences of 'inf' log(2) values by adding 0.001 to each raw value of 0
			$con_val = ($con_val + 0.01) if $con_val == 0; 
			$expt_val = ($expt_val + 0.01) if $expt_val == 0;
			
			my $sign = $expt_val/$con_val;
			
			# Takes the lowest q value for each gene
			# Only include q-values if they are associated with a trascript with expt FPKM >1  
			# as well as having a +ve FC (if considering up-reguulated genes)
			unless ($con eq 'control' ){
				# next unless $q < 0.05;
				$q = 1 if $expt_val < 1 and $direction eq 'up';
				$q = 1 if $sign < 1;
				if ( (not exists $qs{$gene}) or $qs{$gene} > $q ){
				$qs{$gene} = $q;
				}
			}
			
			# Foreach gene/condition sum up all the raw values						
			$con_sum{$gene}{$con} += $con_val;
			$expt_sum{$gene}{$con} += $expt_val;
			
			$con_val = $con_sum{$gene}{$con};
 			$expt_val = $expt_sum{$gene}{$con};	
						
			my $raw_val = $expt_val/$con_val;
	 		$change_sum = eval sprintf('%.2f', (log($raw_val)/ log(2)));
					
			
			# Make a new hash. As the loop cycles though sorted change values, the last value will be the highest count (and thus the correct log(2) value for that gene/condition). Use this in the next loop to pull out only those genes that are all positive or negative...
				
			$genes_by_count{$gene}{$con} = [$xloc, $con_val, $expt_val, $change_sum, $count, $qs{$gene}];
			$all{$gene}{$con} = [$xloc, $con_val, $expt_val, $change_sum, $count, $q];		
								
			# print "$gene\t$con\t$xloc, $con_val, $expt_val, $change_sum, $count\n" if $gene =~ /$query/;
			# push @xlocs, $xloc if $gene =~ /$query/ and $con eq 'six1-eya1';
		
		}
		
		my $num = 2 if $direction eq 'up';
		$num = 1 if $direction eq 'down';
		
		# Only genes with log(2) > threshold allowed through
		# Only genes with chx+dex FPKM > 1 allowed through
		delete $genes_by_count{$gene}{$con} if abs $genes_by_count{$gene}{$con}[3] < $threshold; # remove low FC
		delete $genes_by_count{$gene}{$con} if $genes_by_count{$gene}{$con}[$num] < 1; # remove low FPKM
		delete $genes_by_count{$gene}{$con} if $genes_by_count{$gene}{$con}[5] > 0.05; # remove high q vals
		
		delete $all{$gene}{$con} if $all{$gene}{$con}[$num] < 1;
		$pos_count{$gene}++ if $change_sum > 0 and $con ne 'control';
		$neg_count{$gene}++ if $change_sum < 0 and $con ne 'control';
	}
}

# Experimental conditions only
# Only look for genes that are up regulated in all conditions
for my $gene (keys %genes_by_count) {
	if  (defined $expt_count{$gene} and ( ($pos_count{$gene} or $neg_count{$gene}) == $expt_count{$gene}) ) {
	 	next if exists $neg_count{$gene} and $direction eq 'up';
		next if exists $pos_count{$gene} and $direction eq 'down';
		# Use to filter out specific conditions for 2C lists
		# next if $genes_by_count{$gene}{$filter};

 			foreach my $con (keys $genes_by_count{$gene}) {
				next if $con eq 'control';
 				$seen{$gene}++;
					foreach my $count (  @{$genes_by_count{$gene}{$con}}[4] ){
						my ($xloc, $con_val, $expt_val, $change_val, $count, $q) = @{$genes_by_count{$gene}{$con}};
						# only real duplicates are printed
						$count = '' unless $count > 0;
						$filtered{$gene}{$con} = [$xloc, $con_val, $expt_val, $change_val, $count, $q];
						
						
						# Put the largest change value for each gene in new hash, as long as it's > 1
						# This means that if while all genes log(2) > 0.5 are considered, we need at least 1 of the conditions to have log(2) > 1
						
							if (( (not exists $change{$gene} ) || (abs $change{$gene}[0] < abs $change_val ) ) && (abs $change_val >= $threshold) ) {
								$change{$gene} = [$change_val, $q];
							}
 					}
 			}
 	 }
}


my %control_filt;
for my $gene (keys %all) {
	next unless exists $change{$gene};
  		foreach my $con (keys $all{$gene}) {
 			next unless $con eq 'control';
 				foreach my $count ( @{$all{$gene}{$con}}[4] ) {
 					my ($xloc, $con_val, $expt_val, $change_val, $count, $q) = @{$all{$gene}{$con}};	
 					$count = '' unless $count > 0;
					$control_filt{$gene} = [$con_val, $expt_val, $change_val, $count, $q];												
				}
		}
}	 

my $overlapping_genes = 'overlapping_genes.csv';
open my $overlap_print, '>', $overlapping_genes or die "Can't write to $overlapping_genes: $!";


my $overlapping_count = 0;
my $control_seen = 0;

my %control_ratio;

print "Calculating control ratios for DE genes..\n\n";

fraction(%change);

# read in the transcriptome for two conditions
open my $transcriptome, '<', '/Users/Nick/Documents/Bioinformatics/Data/Merged_reps/transcripts.txt' or die $!;
open my $sequences, '>', '/Users/Nick/Desktop/sequences.txt';

my (%seqs, $xloc);

while (<$transcriptome>){
	chomp;
	($xloc) = />(XLOC_\d+)\./ if /^>/;
 	$seqs{$xloc} = $_ if /[^>]/;
}

# Find Xlocs for queried gene
# open my $xloc_print, '>', "/Users/nickriddiford/Desktop/Analysis/GenesbyXLOC/xlocs.$query.txt" or die $!;
# foreach (@xlocs){
# 		print $xloc_print ">$query($_)\n$sixeya1_seqs{$_}\n" if exists $sixeya1_seqs{$_};
# }

print $overlap_print "Rank,Gene,Control,Six1/Eya1,Fold change log(2),q-val,Control chx,Control chx+dex,Fold change log(2)\n";

my ($x, $control, $eya1mo, $change, $count, $q);

say "Printing ranked DE gene list to 'overlapping_genes.csv'";
say "Printind DE sequences to '/Users/Nick/Desktop/sequences.txt'\n";


# Sort first by rank, and second by log(2)
for my $gene (sort { abs $control_ratio{$b}[1] <=> abs $control_ratio{$a}[1] # FC
				or	 abs $control_ratio{$a}[2] <=> abs $control_ratio{$b}[2] # q-val
				or 	 abs $control_ratio{$a}[0] <=> abs $control_ratio{$b}[0] # rank	  
			 } keys %control_ratio) {
	my ($con_val, $expt_val, $change_val, $count, $q) = @{$control_filt{$gene}} if exists $control_filt{$gene};
	($con_val, $expt_val, $change_val, $count, $q) = (0,0,0, '', 0) if not exists $control_filt{$gene};
	print $overlap_print "@{$control_ratio{$gene}}[0],$gene,";		
		for my $condition (sort { $filtered{$gene}{$a} <=> $filtered{$gene}{$b} } keys %{ $filtered{$gene} } ) { 
			$count++ if $count;
			($x, $control, $eya1mo, $change, $count, $q) = @{$filtered{$gene}{'merge'}} if exists $filtered{$gene}{'merge'};
    	}
		exists $filtered{$gene}{'merge'} ? print $overlap_print "$control,$eya1mo,$change,$q," : print $overlap_print  "-,-,-,";	
		print $overlap_print "$con_val,$expt_val,$change_val\n";
		
		print $sequences "$x $gene\n$seqs{$x}\n";
		# sequence($x, $control, $eya1mo, $change, $gene, $count, 'eya1MO');
		
		# exists $filtered{$gene}{'six1'} ? sequence($Sx, $Sc,$Se,$Sch, $gene, $Sco, 'six1') : sequence($SEx, $SEc,$SEe,$SEch, $gene, $SEco, 'six1-eya1');
}

say "$overlapping_count genes positively upregulated in >= $n conditions";
say "$control_seen genes were removed with Expt.log(2)/Cont.log(2) > 0.5";


my $sleep = 2;
while($sleep--){
    sleep(1);
}

system("open overlapping_genes.csv");


###################
###### SUBS #######
###################

# Potentially quite a uselesss sub. This could be moved to the main sorting loop above.
# sub sequence {
# 	my (%seqs, $xloc);
# 	my ($xlocs, $c, $e, $ch, $gene, $count, $condition) = @_;
# 	$gene =~ s/^ //g;
# 	# $condition eq 'six1' ? print $sequences ">$xlocs $gene $condition \n$six_seqs{$xlocs}\n" : print $sequences ">$xlocs $gene $condition\n$sixeya1_seqs{$xlocs}\n";
# 	print $sequences ">$xlocs $gene $condition \n$seqs{$xlocs}\n"
# }

# 'fraction' takes genes in %change and calculates the ratio between the highest expt. condition and the control value
# and filters out genes with a value > [...]. Outputs %control_ratio to sort genes by this ratio.

sub fraction {
	for my $gene (keys %change){
		next unless $seen{$gene} == $n;
		
		my ($con_val, $expt_val, $change_val, $count, $q) = @{$control_filt{$gene}} if exists $control_filt{$gene};
		($con_val, $expt_val, $change_val, $count) = (0,0,0.001, '') if not exists $control_filt{$gene};
		# Ranks the control genes -ve reg. as higher than absent 
		
		$change_val = 0 if $change_val < 0 and $direction eq 'up'; # For up-reg genes
		$change_val = 0 if $change_val > 0 and $direction eq 'down'; # For down-reg genes
		my $fraction =  eval abs sprintf('%.6f', ($change_val/$change{$gene}[0]));
			
			if  ($fraction > 0.5){
				 $control_seen++;
			 	next;
			}
			else {
				$control_ratio{$gene} = [$fraction, $change{$gene}[0], $change{$gene}[1] ] ;
				$overlapping_count++;
			}	
	}
}

# Extracts information from blast output foreach @ARGV and outputs %gene_list
sub blast_extractor {
    my $file = shift;
    open my $output, '<', $file or die "Can't read file '$file' [$!]\n";
	while (<$output>) { 
		chomp;
		next if /^\s*$/;
		my ($xloc) = /(XLOC_\d+)/;
		my ($change) = /Change:(-?\d+\.\d+|-?inf|-?\d+)/;
		my ($q) = /_q:(\d+\.\d+|-?\d+)/;
		my ($container) = /XLOC_\d+_(.+:.+),_Change/g;
		my ($con_val, $expt_val) = split(/:/, $container);
		my @split = split(/\t/);
		my ($condition, $percent_id, $gene) = ($split[0], $split[2], $split[3]);
		$gene =~ s/_/ /g;
		$gene =~ s/,/ /g;		
		my @small_split = split(/\|/, $gene);
		my $small_gene = $small_split[4];		
		my ($count) = $duplicates{$small_gene}{$condition}++;
		$gene_list{$small_gene}{$condition}{$count} = [$xloc, $con_val, $expt_val, $change, $q];					
	}
	return \%gene_list;
}

sub usage {
    print "\n****co-differenentially expressed genes****\n\n";
    print "Usage: overlapping.pl <de gene list condition 1> .. <de gene list condition N>\n";
    print "This program takes N blast hit results and searches for overlapping genes.\n";
    print "Only conditions with a log(2) fold change value > 0.5 and genes that are upregulated by at least log(2) == 1\n";
    print "are considered.\n\n";
    print "All control values are parsed, and are used to calculate the ratio between control change value/experimental change value.\n";
    print "Output will contain overlapping genes sorted on this ratio, and conditions sorted on change value\n\n";		
}