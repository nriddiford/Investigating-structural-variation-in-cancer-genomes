#!/usr/bin/perl
use warnings;
use strict;

my $file = $ARGV[0];

my ($name, $extention) = split(/\.([^.]+)$/, $file, 2);

open my $in, '<', $file or die $!;

open my $somatic, '>', $name . "_somatic." . $extention;
open my $germline, '>', $name . "_germline." . $extention;
open my $loh, '>', $name . "_LOH." . $extention;

while(<$in>){
    chomp;
	
	if ($_ =~ /^\#/){
		 print $somatic "$_\n";
		 print $germline "$_\n";
		 print $loh "$_\n";
	}
	
	else {
		my ($info) = (split)[7];
		# my ($ss) = (split(/;/, $info))[1];
 		my ($ss) = $info =~ /SS=(\d+)/;
		my ($field) = (split)[4];

		unless ($field =~ /\//){
   			print $germline "$_\n" if $ss == 1;
   			print $somatic "$_\n" if $ss == 2;
   			print $loh "$_\n" if $ss == 3;
		}
	}
}

