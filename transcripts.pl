###################
### Transcripts ###
###################

#!/usr/bin/perl
use warnings;
use strict; 
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

# This script will read in exonic positions from cufflinks output file 'merged.gtf' and relate then to 
# genomic positions. The program then prints out 

# a) All transcript variants for a given gene
# b) All exons in each gene
# c) The longest transcript variant for all assembled genes (numbered by transcript variant)

# Run as perl transcripts.pl in the same directory as 'merged.gtf'

open my $merged, '<', 'merged.gtf' or die "usage perl transcripts.pl\nIs 'merged.gtf'in the current directory?\n";

my (%gene, $scaf, $start, $stop, $xloc, $tcons, $exon, $length, $info);
 
while(<$merged>){
	my @split = split('\t');
	($scaf) = ">$split[0]";
	($start) = $split[3];
	($stop) = $split[4];
	($info) = $split[8];
	($xloc, $tcons, $exon) = $info =~ /gene_id "(.+?)"; transcript_id "(.+?)"; exon_number "(\d+)";/g;
	$gene{$xloc}{$tcons}{$exon} = [$start, $stop, $scaf];
}

open my $genome, '<', '/Users/Nick/Documents/Bioinformatics/Genomes/Laevis_7/Xl7.fa' or die $!;

my (%seq, $head);
while (<$genome>) {
  chomp;
  $head = $_ if /^>/;   
  $seq{$head} = $_ if /^[A-Z]/;
}

open my $transcript_variant, '>', 'transcript_variants.txt' or die $!;
open my $exons, '>', 'transcripts_by_exon.txt' or die $!;
open my $transcripts, '>', 'transcripts.txt' or die $!;

my ($exon_seq, %transcript, %seq_length, $exon_length);

my %most_exons;
my $splice_var;

for my $xloc (sort keys %gene){
	my $seen = 0;
	for my $tcons (sort keys $gene{$xloc}){
		$seen++;
		$splice_var = "$xloc.$seen";
		my $exon_count = 0;
			for my $exon (sort { $a <=> $b } keys $gene{$xloc}{$tcons}){
				my ($start, $stop, $scaf) = @{$gene{$xloc}{$tcons}{$exon}};
				$exon_count++;
				$start = ($start - 1);
				$exon_length = $stop - $start;
				$exon_seq = substr $seq{$scaf}, $start, $exon_length;
				print $exons ">$xloc.$seen\_$exon\n$exon_seq\n";	
				$transcript{$tcons} .= $exon_seq;					
					if ( (not exists $most_exons{$xloc}) or ($most_exons{$xloc}[0] < $exon_count) ){
						$most_exons{$xloc} = [$exon_count, $splice_var, $transcript{$tcons}];
					}			
		}
		print $transcript_variant ">$xloc.$seen\n$transcript{$tcons}\n";	
	}	
}

for my $gene (sort keys %most_exons){
	my ($exon_count, $splice_var, $sequence) = @{$most_exons{$gene}};
	print $transcripts ">$splice_var\n$sequence\n";
}



