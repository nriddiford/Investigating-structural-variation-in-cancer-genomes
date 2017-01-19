#!/usr/bin/perl
use warnings;
use strict;

use Data::Dumper;
use feature qw/ say /;

use Bio::SeqIO;

my $file = $ARGV[0];

exit usage() unless $#ARGV == 0;

my $seqio = Bio::SeqIO->new(-file => $file,
							-format => 'Fasta'
							);
							
# my $seq_out = Bio::SeqIO->new( -file  => "genome.index.fa",
# 							   -format => 'Fasta',
#                             );
						
open my $out, '>', 'genome.index.txt' or die $!;

# my ($name, $extention) = split(/\.([^.]+)$/, $file, 2);

my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

# my %chroms = qw / X 23542271 /;

while(my $seq = $seqio->next_seq) {
	my $header = $seq->primary_id;
	my $string = $seq->seq;
	
	if ($chroms{$header}){
	
		for my $i ( 1 .. length $string ){
			# $i -1 to get 1-based co-ordinates
	 		my ($nuc) = substr( $string, $i - 1, 1 );
	 		print $out "$header\t$nuc\t$i\n";
		}
	
	}

}

sub usage {
	say "Usage: perl $0 <genome.fa>";
}

