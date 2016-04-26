################
### Godzilla ###
################

#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use feature qw(say);
use autodie;

# Run `perl godzilla.pl` for usage statemnent
# This program will take N blast hit results and output those genes which are found to be overlapping in $n conditions.
# Run as follows: `perl godzilla.pl control.txt six.txt six-eya.txt eya.txt` where blast output has the following, tab delimited format:

# gene	XLOC_control.fpkm:expt.fpkm,_Change:FC_val	%id	blast hit

# Always run 'control' condition first

if ($#ARGV < 0) {
	usage();
	exit;
}

# set a flag for 'silent run' which will not read in transcriptomes or open .csv file behind script
my $flag1;
if ($ARGV[0] eq '-s'){
	$flag1 = shift @ARGV;
}
else {$flag1 = 1;}

my (%duplicates, %gene_list, %change, %seen, %seen_all);

print "Opening blast ouput foreach condition...\n";

blast_extractor($_) foreach @ARGV;# To extract data from blast output

# Select for up or down regulation
print "Enter direction: ";
chomp(my $direction = <STDIN>);

# Select the condition being run to read in corresponding transcriptome
my $seqs_to_print;
print "Enter sequences to print (six1 = 0 eya1 = 1 six1-eya1 = 2): ";
chomp($seqs_to_print = <STDIN>);

my (%filtered, %expt_count, %pos_count, %neg_count, %con_sum, %expt_sum, %genes_by_count, %all, @xlocs);
my $change_sum;

my $threshold = 0.5;
my $n = 1;

say "Filtering genes for DE genes...";

my $gene;
foreach $gene (sort keys %gene_list) {
	foreach my $con (sort keys $gene_list{$gene}) {
		$expt_count{$gene}++ unless $con eq 'control';
		foreach my $count (sort keys $gene_list{$gene}{$con} ){
			my ($xloc, $con_val, $expt_val, $change_val) = @{$gene_list{$gene}{$con}{$count}};
			# Remove the occurrences of 'inf' log(2) values by adding 0.001 to each raw value of 0
			$con_val = ($con_val + 0.001) if $con_val == 0;
			$expt_val = ($expt_val + 0.001) if $expt_val == 0;

			# Foreach gene/condition sum up all the raw values
			$con_sum{$gene}{$con} += $con_val;
			$expt_sum{$gene}{$con} += $expt_val;

			$con_val = $con_sum{$gene}{$con};
 			$expt_val = $expt_sum{$gene}{$con};

			my $raw_val = $expt_val/$con_val;
	 		$change_sum = eval sprintf('%.6f', (log($raw_val)/ log(2)));

			# Make a new hash. As the loop cycles though sorted change values, the last value will be the highest count (and thus the correct log(2) value for that gene/condition). Use this in the next loop to pull out only those genes that are all positive or negative...

			$genes_by_count{$gene}{$con} = [$xloc, $con_val, $expt_val, $change_sum, $count];
			$all{$gene}{$con} = [$xloc, $con_val, $expt_val, $change_sum, $count];

		}
		# Only genes with log(2) > 0.5 allowed through
		# Only genes with chx+dex FPKM > 1 allowed through
		
		# In upeg run mode, consider genes with chx+dex FPKM > 1 as expressed 
		# In downreg run mode, consider genes with chx FPKM > 1 as expressed
		
		my $num = 2 if $direction eq 'up';
		$num = 1 if $direction eq 'down';

		delete $genes_by_count{$gene}{$con} if abs $change_sum < $threshold or $genes_by_count{$gene}{$con}[$num] < 1;
		# Treat control in the same way
		delete $all{$gene}{$con} if $all{$gene}{$con}[$num] < 1 and $con eq 'control';

		$pos_count{$gene}++ if $change_sum > 0 and $con ne 'control';
		$neg_count{$gene}++ if $change_sum < 0 and $con ne 'control';
	}
}

# Experimental conditions only
# Only look for genes that are up or down reg in all conditions
for my $gene (sort keys %genes_by_count) {
	if  ( defined $expt_count{$gene} and ( ($pos_count{$gene} or $neg_count{$gene}) == $expt_count{$gene}) ) {
		next if exists $neg_count{$gene} and $direction eq 'up';
		next if exists $pos_count{$gene} and $direction eq 'down';
 			foreach my $con ( sort keys $genes_by_count{$gene}) {
				next if $con eq 'control';
 				$seen{$gene}++;
					foreach my $count (  @{$genes_by_count{$gene}{$con}}[4] ){
						my ($xloc, $con_val, $expt_val, $change_val, $count) = @{$genes_by_count{$gene}{$con}};
						# only real duplicates are printed
						$count = '' unless $count > 0;
						$filtered{$gene}{$con} = [$xloc, $con_val, $expt_val, $change_val, $count];

						# Put the smallest change value for each gene in new hash, as long as it's > 1
						# This means that if while all genes log(2) > 0.5 are considered, we need at least 1 of the conditions to have log(2) > 1

							if (( (not exists $change{$gene} ) || (abs $change{$gene} < abs $change_val ) ) && (abs $change_val >= $threshold) ) {
								$change{$gene} = $change_val;
							}
 					}
 			}
 	 }
}

my %control_filt;
for my $gene (sort keys %all) {
	next unless exists $change{$gene};
		foreach my $con (sort keys $all{$gene}) {
 			next unless $con eq 'control';
 				foreach my $count ( @{$all{$gene}{$con}}[4] ) {
 					my ($xloc, $con_val, $expt_val, $change_val, $count) = @{$all{$gene}{$con}};
 					$count = '' unless $count > 0;
					$control_filt{$gene} = [$con_val, $expt_val, $change_val, $count];
				}
		}
}

my $overlapping_count = 0;
my $control_seen = 0;
my %control_ratio;

print "Calculating control ratios for DE genes..\n\n";

fraction(%change);

my $sequences;
my (%six_seqs, %eya_seqs, %sixeya1_seqs, $xloc);

# Skip reading in transcriptomes if in silent mode
unless ($flag1 eq '-s'){

open $sequences, '>', "/Users/nickriddiford/Desktop/$seqs_to_print.$direction.sequences.txt" or die $!;

if ($seqs_to_print == 0){
	print "Opening transcriptome for six1...\n";
	open my $S_transcriptome, '<', '/Users/nickriddiford/Desktop/Data/Genes/six1/transcripts.txt' or die $!;
	while (<$S_transcriptome>){
		chomp;
		($xloc) = />(XLOC_\d+)\./ if /^>/;
 		$six_seqs{$xloc} = $_ if  /[^>]/;
	}
}
if ($seqs_to_print == 1){
	print "Opening transcriptome for eya1...\n";
	open my $E_transcriptome, '<', '/Users/nickriddiford/Desktop/Data/Genes/eya1/transcripts.txt' or die $!;
	while (<$E_transcriptome>){
		chomp;
		($xloc) = />(XLOC_\d+)\./ if /^>/;
 		$eya_seqs{$xloc} = $_ if  /[^>]/;
	}
}
if ($seqs_to_print == 2){
	print "Opening transcriptome for six1-eya1...\n";
	open my $SE_transcriptome, '<', '/Users/nickriddiford/Desktop/Data/Genes/six1-eya1/transcripts.txt' or die $!;
	while (<$SE_transcriptome>){
		chomp;
		($xloc) = />(XLOC_\d+)\./ if /^>/;
 		$sixeya1_seqs{$xloc} = $_ if /[^>]/;
	}
}

}

my $overlapping_genes = 'overlapping_genes.csv';
open my $overlap_print, '>', $overlapping_genes or die "Can't write to $overlapping_genes: $!";

print $overlap_print "Rank,Gene,Six1 CHX,Six1 CHX+DEX,Log(2),Eya1 CHX,Eya1 CHX+DEX,Log(2),Six1-Eya1 CHX,Six1-Eya1 CHX+DEX,Log(2),Control CHX,Control CHX+DEX,Log(2)\n";

my ($Sx, $Sc, $Se, $Sch, $Sco, $Ex, $Ec, $Ee, $Ech, $Eco, $SEx, $SEc, $SEe, $SEch, $SEco);

say "Printing ranked DE gene list to 'overlapping_genes.csv'";
say "Printind DE sequences to '/Users/nickriddiford/Desktop/sequences.txt'\n";

# Output 'out.txt' as input for 'sets.pl'
open my $out, '>', "/Users/nickriddiford/Desktop/$seqs_to_print.$direction.out.txt";

# Sort first by rank, and second by FC
for my $gene (sort { abs $control_ratio{$a}[0] <=> abs $control_ratio{$b}[0] 
					or
					abs $control_ratio{$b}[1] <=> abs $control_ratio{$a}[1] } 
					keys %control_ratio) {
	my ($con_val, $expt_val, $change_val, $count) = @{$control_filt{$gene}} if exists $control_filt{$gene};
	($con_val, $expt_val, $change_val, $count) = (0,0,0, '') if not exists $control_filt{$gene};
	print $overlap_print "$control_ratio{$gene}[0],$gene,";

	print $out "$control_ratio{$gene}[0]\t$gene\t";
	
		for my $condition (sort { $filtered{$gene}{$a} <=> $filtered{$gene}{$b} } keys %{ $filtered{$gene} } ) {
			$count++ if $count;
			($Sx, $Sc, $Se, $Sch, $Sco) = @{$filtered{$gene}{'six1'}} if exists $filtered{$gene}{'six1'};
			($Ex, $Ec, $Ee, $Ech, $Eco) = @{$filtered{$gene}{'eya1'}} if exists $filtered{$gene}{'eya1'};
			($SEx, $SEc, $SEe, $SEch, $SEco) = @{$filtered{$gene}{'six1-eya1'}} if exists $filtered{$gene}{'six1-eya1'};
    	}
		# Input for 'sets.pl'
		print $out "$Sx\t" if exists $filtered{$gene}{'six1'};
		print $out "$Ex\t" if exists $filtered{$gene}{'eya1'};
		print $out "$SEx\t" if exists $filtered{$gene}{'six1-eya1'};
		
		exists $filtered{$gene}{'six1'} ? print $overlap_print "$Sc,$Se,$Sch," : print $overlap_print  "-,-,-,";
		exists $filtered{$gene}{'eya1'} ? print $overlap_print "$Ec,$Ee,$Ech," : print $overlap_print  "-,-,-,";
		exists $filtered{$gene}{'six1-eya1'} ? print $overlap_print "$SEc,$SEe,$SEch," : print $overlap_print  "-,-,-,";
		print $overlap_print "$con_val,$expt_val,$change_val\n";
		
		# Input for 'sets.pl'
		exists $filtered{$gene}{'six1'} ? print $out "$Sc\t$Se\t$Sch\t" : print $out  "-\t-\t-\t";
		exists $filtered{$gene}{'eya1'} ? print $out "$Ec\t$Ee\t$Ech\t" : print $out  "-\t-\t-\t";
		exists $filtered{$gene}{'six1-eya1'} ? print $out "$SEc\t$SEe\t$SEch\t" : print $out  "-\t-\t-\t";
		print $out "$con_val\t$expt_val\t$change_val\n";
		
		unless ($flag1 eq '-s'){
			sequence($Sx, $Sc, $Se, $Sch, $gene, $Sco, 'six1') if $seqs_to_print == 0;
	 		sequence($Ex, $Ec, $Ee, $Ech, $gene, $Eco, 'eya1') if $seqs_to_print == 1;
			sequence($SEx, $SEc, $SEe, $SEch, $gene, $SEco, 'six1-eya1')if $seqs_to_print == 2;
		}
}

say "$overlapping_count genes $direction-regulated in $n conditions";
say "$control_seen genes were removed with Expt.log(2)/Cont.log(2) > 0.5";


unless ($flag1 eq '-s'){
my $sleep = 2;
while($sleep--){
    sleep(1);
}

system("open overlapping_genes.csv");
}
###################
###### SUBS #######
###################

# Attach sequences
unless ($flag1 eq '-s'){
	sub sequence {
		my (%seqs, $xloc);
		my ($xlocs, $c, $e, $ch, $gene, $count, $condition) = @_;
		$gene =~ s/^ //g;
		print $sequences ">$xlocs $gene $condition \n$six_seqs{$xlocs}\n" if $seqs_to_print == 0;
		print $sequences ">$xlocs $gene $condition \n$eya_seqs{$xlocs}\n" if $seqs_to_print == 1;
		print $sequences ">$xlocs $gene $condition \n$sixeya1_seqs{$xlocs}\n" if $seqs_to_print == 2;
	}
}

# 'fraction' takes genes in %change and calculates the ratio between the highest expt. condition and the control value
# and filters out genes with a value > 0.5.
# Outputs %control_ratio to sort genes by this ratio.

sub fraction {
	for my $gene ( sort keys %change){
		next unless $seen{$gene} == $n;

		my ($con_val, $expt_val, $change_val, $count) = @{$control_filt{$gene}} if exists $control_filt{$gene};
		($con_val, $expt_val, $change_val, $count) = (0,0,0.001, '') if not exists $control_filt{$gene};
		
		# Ranks the control genes -ve reg. as higher than absent
		if ($direction eq 'up'){
			$change_val = 0 if $change_val < 0;
		}
		else {$change_val = 0 if $change_val > 0;}

		my $fraction =  eval abs sprintf('%.3f', ($change_val/$change{$gene}));

			if (eval sprintf('%.1f', $fraction) > 0.5){
				 $control_seen++;
			 	next;
			}
			else {
				$control_ratio{$gene} = [$fraction, $change{$gene}] ;
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
		my ($container) = /XLOC_\d+_(.+:.+),_Change/g;
		my ($con_val, $expt_val) = split(/:/, $container);
		my @split = split(/\t/);
		my ($condition, $percent_id, $gene) = ($split[0], $split[2], $split[3]);
		$gene =~ s/_/ /g;
		$gene =~ s/,/ /g;
		my @small_split = split(/\|/, $gene);
		my $small_gene = $small_split[4];
		my ($count) = $duplicates{$small_gene}{$condition}++;
		$gene_list{$small_gene}{$condition}{$count} = [$xloc, $con_val, $expt_val, $change];
	}
	return \%gene_list;
}

sub usage {
    print "\n****GODZILLA****\n\n";
    print "Usage: godzilla.pl <control gene list> <de gene list condition 1> .. <de gene list condition N>\n";
    print "This program takes N blast hit results and searches for overlapping genes.\n";
    print "Only conditions with a log(2) fold change value > $threshold are considered.\n\n";
    print "All control values are parsed, and are used to calculate the ratio between control change value/experimental change value.\n";
    print "Output in .csv will contain genes sorted first on rank and second on FC\n\n";
}