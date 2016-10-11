#################
### annotator ###
#################

#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;


# This script takes blast results as input and annotates all transcripts with best hit and user specified condition name
# It also counts the number of genes that return no hit. Output as two .txt files which can be blasted against a different database
# (or with a lower E-value) and added to the blast output.

# Run as perl annotator.pl <Xenopus blast output> <condition> <de_genes.txt>

if ($#ARGV != 2){
	usage();
	exit
}

my $input1 = $ARGV[0];
open my $blast, '<', $input1 or die $!;

my $condition = $ARGV[1];

# edit.input will trim input, removing genes that don't hit the database
my $editedblast = "$condition.$input1";
open my $unknowns_removed, '>', $editedblast or die "Can't write to $editedblast: $!";

# Split output from blast and outout edited input (with non-hits removed)
my $blast_count = 0;
my (%blasthits, @test, @blastID);
while (<$blast>) {
	chomp;
	$blast_count++;
	my @genematch = ( $_ =~ /(XLOC_\d+)/g);
	push @blastID, @genematch;
	my @split = split('\s+');
	print $unknowns_removed "$condition\t$_\n";

}

push @{$blasthits{$blastID[$_] } }, [ $blastID[$_] ] for 0 .. $#blastID;

# read in de_genes
my $de_seqs = $ARGV[2];
open my $de_genes, '<', $de_seqs or die "Can't open '$de_seqs'";

my $gene_count = 0;
my (@gene, @head, @seq);
while (<$de_genes>) {
	chomp;
	my @match = ( $_ =~ /(XLOC_\d+)/g);
	push @gene, @match;
	push @head, $_ if /^>/;
	push @seq, $_ if /^[A-Z]/;
	$gene_count++ if /^>/;
}

print "$gene_count queries\n";
print "$blast_count hits in $input1 \n";

my %genes;
push @{ $genes{ $gene[$_] } }, [ $head[$_],  $seq[$_] ] for 0 .. $#gene;

my $nohit = 'unknown_genes.txt';
open my $unknown_print, '>', $nohit or die "Can't write to $nohit: $!";

# Output to 'unknown_genes.txt' all the de_genes that return no blast hit
my $i = 0;
for my $key (sort keys %genes) {
	unless (exists $blasthits{$key}) {
		$i++;
		for my $part (@ { $genes{$key} } ) {
				my ($head, $seq) = @$part;
				print $unknown_print "$head\n$seq\n";
		}
	}
}
print "$i genes that return no blast hit\n";

##########################################

sub usage {
    print "\n****ANNOTATOR****\n\n";
    print "Usage: annotator.pl <Xenopus blast output> <condition> <de_genes.txt>\n\n";
    print "Takes BLAST results as input and finds all genes that return no hit (not present\n";
    print "in the blast output but present in the de_genes).\n\n";
    print "Outputs: 'unknown_genes.txt' to BLAST against a different database\n";
    print "Outputs: 'edit.<input>.txt' to be used as input for 'godzilla.pl'\n\n";
    print "Nick Riddiford\n";
    print "August 2013\n\n";
}
