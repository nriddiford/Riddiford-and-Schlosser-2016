############
### Sets ###
############

#!/usr/bin/perl
use warnings;
use strict; 
use Spreadsheet::ParseXLSX;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

# This script divides blast input into discrete sets. 
# The order input is given disctates the naming of the output: 
# if run as perl sets.pl list1 list2 list3
# output will be as follows:

# SetA = list1 not list2 not list3
# SetB = list1 and list2 not list3
# SetC = list1 and list3 not list2
# SetD = list1 and list2 and list3

unless ($#ARGV == 2) {
	usage();
	exit;
}

# used to specify which sequence to print (and which column to sort on)
print "Enter primary set: ";
chomp(my $set = <STDIN>);

my (%c1_vals, %c2_vals, %c3_vals);

parse($_) foreach @ARGV;

my $count = 0;
my $c;
my ($c1, $c2, $c3);

my (%c1_data, %c2_data, %c3_data);


 
my $c1_c2_count = 0;
my $c1_c3_count = 0;
my $cond1 = 0;
my $total = 0;
my $all = 0;

my (%six_seqs, %eya_seqs, %sixeya1_seqs, $xloc);

open my $S_transcriptome, '<', '/Users/nickriddiford/Desktop/Data/Genes/six1/transcripts.txt' or die $!;
while (<$S_transcriptome>){
	chomp;
	($xloc) = />(XLOC_\d+)\./ if /^>/;
 	$six_seqs{$xloc} = $_ if  /[^>]/;
}

open my $E_transcriptome, '<', '/Users/nickriddiford/Desktop/Data/Genes/eya1/transcripts.txt' or die $!;
while (<$E_transcriptome>){
	chomp;
	($xloc) = />(XLOC_\d+)\./ if /^>/;
	$eya_seqs{$xloc} = $_ if  /[^>]/;
}

open my $SE_transcriptome, '<', '/Users/nickriddiford/Desktop/Data/Genes/six1-eya1/transcripts.txt' or die $!;
while (<$SE_transcriptome>){
	chomp;
	($xloc) = />(XLOC_\d+)\./ if /^>/;
	$sixeya1_seqs{$xloc} = $_ if /[^>]/;
}

open my $set_A_seqs, '>', "$c1\_set_A_seqs.txt" or die $!;
open my $set_B_seqs, '>', "$c1\_set_B_seqs.txt" or die $!;
open my $set_C_seqs, '>', "$c1\_set_C_seqs.txt" or die $!;
open my $set_D_seqs, '>', "$c1\_set_D_seqs.txt" or die $!;


my $switch;
$switch = 5 if $set eq 'six';
$switch = 8 if $set eq 'eya';
$switch = 11 if $set eq 'six-eya';

open my $space_A, '>', "$c1\_set_A.csv" or die $!;
open my $space_B, '>', "$c1\_set_B.csv" or die $!;
open my $space_C, '>', "$c1\_set_C.csv" or die $!;
open my $space_D, '>', "$c1\_set_D.csv" or die $!;

print $space_A "Rank,Gene,Six1 CHX,Six1 CHX+DEX,Log(2),Eya1 CHX,Eya1 CHX+DEX,Log(2),Six1-Eya1 CHX,Six1-Eya1 CHX+DEX,Log(2),Control CHX,Control CHX+DEX,Log(2)\n"; 
print $space_B "Rank,Gene,Six1 CHX,Six1 CHX+DEX,Log(2),Eya1 CHX,Eya1 CHX+DEX,Log(2),Six1-Eya1 CHX,Six1-Eya1 CHX+DEX,Log(2),Control CHX,Control CHX+DEX,Log(2)\n"; 
print $space_C "Rank,Gene,Six1 CHX,Six1 CHX+DEX,Log(2),Eya1 CHX,Eya1 CHX+DEX,Log(2),Six1-Eya1 CHX,Six1-Eya1 CHX+DEX,Log(2),Control CHX,Control CHX+DEX,Log(2)\n"; 
print $space_D "Rank,Gene,Six1 CHX,Six1 CHX+DEX,Log(2),Eya1 CHX,Eya1 CHX+DEX,Log(2),Six1-Eya1 CHX,Six1-Eya1 CHX+DEX,Log(2),Control CHX,Control CHX+DEX,Log(2)\n"; 

# sort the lists first by control rank, then by FC and print out as in 'overlapping.pl'
for my $gene (sort { abs @{$c1_data{$a}}[0]  <=> abs @{$c1_data{$b}}[0]
				or 	abs @{$c1_data{$b}}[$switch] <=> abs @{$c1_data{$a}}[$switch]
					} keys %c1_data) {
						# second sort:  six1 = 5
						#				eya1 = 8
						#				six-eya = 11
	my ($rank, $gene, $xloc, $sc, $scd, $sch, $ec, $ecd, $ech, $sec, $secd, $sech, $cc, $ccd, $cch) = @{$c1_data{$gene}};
	$total++;
	if (not $c2_vals{$gene} and not $c3_vals{$gene}){
		$cond1++;
		print $space_A "$rank,$gene,$sc,$scd,$sch,$ec,$ecd,$ech,$sec,$secd,$sech,$cc,$ccd,$cch\n";
		print $set_A_seqs ">$xloc\_$gene\_$c1\n$six_seqs{$xloc}\n" if $set eq 'six';
		print $set_A_seqs ">$xloc\_$gene\_$c1\n$eya_seqs{$xloc}\n" if $set eq 'eya';
		print $set_A_seqs ">$xloc\_$gene\_$c1\n$sixeya1_seqs{$xloc}\n" if $set eq 'six-eya';
		
	}
	if ($c2_vals{$gene} and not $c3_vals{$gene}){
		$c1_c2_count++;
		print $space_B "$rank,$gene,$sc,$scd,$sch,$ec,$ecd,$ech,$sec,$secd,$sech,$cc,$ccd,$cch\n";
		print $set_B_seqs ">$xloc\_$gene\_$c1\n$six_seqs{$xloc}\n" if $set eq 'six';
		print $set_B_seqs ">$xloc\_$gene\_$c1\n$eya_seqs{$xloc}\n" if $set eq 'eya';
		print $set_B_seqs ">$xloc\_$gene\_$c1\n$sixeya1_seqs{$xloc}\n" if $set eq 'six-eya';
	}
	if ($c3_vals{$gene} and not $c2_vals{$gene}){
		$c1_c3_count++;
		print $space_C "$rank,$gene,$sc,$scd,$sch,$ec,$ecd,$ech,$sec,$secd,$sech,$cc,$ccd,$cch\n";
		print $set_C_seqs ">$xloc\_$gene\_$c1\n$six_seqs{$xloc}\n" if $set eq 'six';
		print $set_C_seqs ">$xloc\_$gene\_$c1\n$eya_seqs{$xloc}\n" if $set eq 'eya';
		print $set_C_seqs ">$xloc\_$gene\_$c1\n$sixeya1_seqs{$xloc}\n" if $set eq 'six-eya';
	}
	if ($c2_vals{$gene} and $c3_vals{$gene}){
		$all++;
		print $space_D "$rank,$gene,$sc,$scd,$sch,$ec,$ecd,$ech,$sec,$secd,$sech,$cc,$ccd,$cch\n";
		print $set_D_seqs ">$xloc\_$gene\_$c1\n$six_seqs{$xloc}\n" if $set eq 'six';
		print $set_D_seqs ">$xloc\_$gene\_$c1\n$eya_seqs{$xloc}\n" if $set eq 'eya';
		print $set_D_seqs ">$xloc\_$gene\_$c1\n$sixeya1_seqs{$xloc}\n" if $set eq 'six-eya';
	}
}
# print numbers for each set 
print "$c1 total = $total\n";
print "$c1 = $cond1\n";
print "$c1 + $c2  = $c1_c2_count\n";
print "$c1 + $c3 = $c1_c3_count\n";
print "$c1+$c2+$c3 = $all\n";
my $check = ($cond1 + $c1_c2_count + $c1_c3_count + $all);
print "check = $check\n";

# read in each file (output from 'godzilla.pl'). $ARGV[0] is set as the 'primary' set (ie that which intersecting lists are found for)

sub parse {
	$count++;
    my $file = shift;
	($c = $file) =~ s/\.[^.]+$//;
    open my $list, '<', $file or die "Can't read file '$file' [$!]\n";
	while(<$list>) {
		chomp;
		if ($count == 1){
			my @split = split(/\t/);
			$c1_vals{$split[1]}++;
			$c1_data{$split[1]} = [@split];
			 $c1 = $c;
		 }
 		if ($count == 2){
			my @split = split(/\t/);
 			$c2_vals{$split[1]}++;
			$c2_data{$split[1]} = [@split];
 			 $c2 = $c;
 		 }
 		if ($count == 3){
			my @split = split(/\t/);
 			$c3_vals{$split[1]}++;
			$c3_data{$split[1]} = [@split];
 			 $c3 = $c;
 		 }
	}
}

sub usage {
	print "Usage: set.pl <gene list1> <gene list2> <gene list3>\n";
	print "Calculates intersections between different sets\n";
}
