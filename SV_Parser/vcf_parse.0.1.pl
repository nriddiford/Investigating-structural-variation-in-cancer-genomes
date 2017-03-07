#!/urs/bin/perl

use VCF_1_0;
use strict;
use warnings;

use feature qw/ say /;
use Data::Dumper;

use Getopt::Long qw/ GetOptions /;

my $vcf_file; 
my $help;
my $id;
my $dump;
my $filter;
my $chromosome;

# Should add score threshold option
GetOptions( 'vcf=s'	        	=>		\$vcf_file,
			'id=s'				=>		\$id,
			'dump'				=>		\$dump,
			'filter'			=>		\$filter,
			'chromosome=s'		=>		\$chromosome,
			'help'              =>      \$help
	  ) or die usage();

if ($help) { exit usage() } 

if (not $vcf_file) {
	 exit usage();
} 


# Retun SV and info hashes 
my ( $SVs, $info ) = VCF_1_0::typer($vcf_file);

# Print all infor for specified id


VCF_1_0::summarise_variants( $SVs, $filter ) unless $id or $dump;

# Print all infor for specified id

VCF_1_0::get_variant( $id, $SVs, $info, $filter ) if $id;

# Dump all variants to screen
VCF_1_0::dump_variants( $SVs, $info, $filter, $chromosome ) if $dump;


sub usage {
	say "********** VCF_parser ***********";
    say "Usage: $0 [options]";
	say "  --vcf = VCF file for parsing";
	say "  --id = extract information for a given variant";
	say "  --dump = cycle through all variants (can be combined with both -f and -c)";
	say "  --filter = apply filters and mark filtered variants";
	say "  --chromosome = used in conjunction with --dump will cycle though variants on chromosome speciified in -c";
	say "  --help\n";
	say "Nick Riddiford 2017";
}