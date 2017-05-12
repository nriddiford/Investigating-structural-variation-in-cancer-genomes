#!/usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Data::Dumper;

exit usage() unless $#ARGV == 2;

my $file = $ARGV[0];
my $FC = $ARGV[1];
my $condition = $ARGV[2];

my ($name, $extention) = split(/\.([^.]+)$/, $file, 2);

open my $in, '<', $file or die $!;

open my $gtf_up, '>', $name . "_topchange_up_" . $FC . ".gff3";
open my $gtf_down, '>', $name . "_topchange_down_" . $FC . ".gff3";
open my $out, '>', $name . "_topchange_" . $FC . ".cnv";

# Print headers to both files
print $gtf_up "##gff-version 3\n";
print $gtf_up "#track name=\"$condition HIGH\" gffTags=on\n";

print $gtf_down "##gff-version 3\n";
print $gtf_down "#track name=\"$condition LOW\" gffTags=on\n";

while(<$in>){
        chomp;
        if ($. == 1){
                print $out "$_\n";
                next;
        }

        #log2 is in field[6]
        my ($val) = (split)[6];
        next if $val eq 'NA';
        # Save absolute (ie both 1 and -1 evaluate as 1)
        my $abs_log2 = abs($val);

        next if $val =~ /-?Inf/;

        if ( $abs_log2 >= $FC ){

            print $out "$_\n";
            my $col;

            my ($chr, $start, $stop, $p, $size) = (split)[0,1,2,7,9];
            $chr =~ s/\"//g;
            $val = sprintf('%.1f', $val);

            if ($val < 0){

		$col = "#608000";
                $col = "#608000" if $abs_log2 >= 1 and $abs_log2 < 2;
                $col = "#808000" if $abs_log2 >= 2 and $abs_log2 < 3;
                $col = "#804000" if $abs_log2 >= 3 and $abs_log2 < 5;
                $col = "#800000" if $abs_log2 >= 5;

                print $gtf_down join ("\t", $chr, ".", "region", $start, $stop, ".", "+", ".", "Name=$val;colour=$col;");
                print $gtf_down "\n";
            }

            elsif ($val > 0){
		
		$col = "#008080";
                $col = "#008080" if $abs_log2 >= 1 and $abs_log2 < 2;
                $col = "#004080" if $abs_log2 >= 2 and $abs_log2 < 3;
                $col = "#000080" if $abs_log2 >= 3 and $abs_log2 < 5;
                $col = "#400080" if $abs_log2 >= 5;

                print $gtf_up join ("\t", $chr, ".", "region", $start, $stop, ".", "+", ".", "Name=$val;colour=$col;");
                print $gtf_up "\n";
            }

         }
}





sub usage {
	say "Usage: $0 <.cnv file> <FC threshold> <gff trackname>";
}
