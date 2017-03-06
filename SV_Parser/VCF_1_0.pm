#!/urs/bin/perl

package VCF_1_0;
use strict;
use warnings;

use feature qw/ say /;
use Data::Dumper;

sub typer {
	my $file = shift;
	my $type;
	if (`grep "source=LUMPY" $file`){
		say "Recognised $file as Lumpy input";
		$type = 'lumpy';
		parse($file, $type);
	}
	
	elsif (`grep "DELLY" $file`){
		say "Recognised $file as Delly input";
		$type = 'delly';
		parse($file, $type);
	}
	
	else {
		die "This VCF is not from lumpy or delly. Abort";
	}
}

sub parse {
	my ($file, $type) = @_;
	open my $in, '<', $file or die $!;

	my @header;
	my %SVs;
	my ($tumour_name, $control_name);
	my %info;
	my (@format_long, @info_long, @info_available);

	my %info_block;
	my $filter_count;

	my %sample_info;
	my @samples;
	
	while(<$in>){
		chomp;

		if (/^#{2}/){
	 		push @header, $_;

			if (/##FORMAT/){

				push @format_long, $_ =~ /\"(.*?)\"/;
			}

			if (/##INFO/) {
				my ($info_long) = $_ =~ /\"(.*?)\"/;
				my ($available_info) = $_ =~ /ID=(.*?),/;
				$info_block{$available_info} = $info_long;
			}
			next;

		}

		if (/^#{1}/){
			my @split = split;
			push @samples, $_ foreach @split[9..$#split];
 			$tumour_name 	= (split)[9];
 			$control_name 	= (split)[10];
			next;
		}

		my @fields = split;
	    my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block) = @fields;


	    my @tumour_parts 	= split(/:/, $tumour_info_block);
	    my @normal_parts 	= split(/:/, $normal_info_block);
		my @format 		 	= split(/:/, $format_block);
		
		my @info_parts		= split(/;/, $info_block);
		
		foreach(@info_parts){
			my ($info_key, $info_value);

			if (/=/){
				($info_key, $info_value) = $_ =~ /(.*)=(.*)/;
			}

			else {
				($info_key) = $_ =~ /(.*)/;
				$info_value = "TRUE";
			}
			$sample_info{$info_key} = $info_value;
		}
		
	    my ($SV_type) = $info_block =~ /SVTYPE=(.*?);/;
	    my ($stop) = $info_block =~ /;END=(.*?);/;
		
		my $non_hom_control_flag = 0;
		
		
		my ($t_GT, $c_GT);
		
		my @filter_reasons;
		
		for (my $i = 0; $i <= $#format; $i++ ){
			if ( $format[$i] eq 'GT' ){
				$t_GT = $tumour_parts[$i];
								
				if ( $normal_parts[$i] eq '1/1' or $normal_parts[$i] eq '0/1' ){
					$non_hom_control_flag = 1;
					push @filter_reasons, 'control_gt=' . $normal_parts[$i];
									
				}
			}
		}
		
		
							
		my ($SV_length, $chr2, $t_SR, $t_PE);
		if ($type eq 'lumpy'){
			( $SV_length, $chr2, $t_SR, $t_PE, @filter_reasons ) = lumpy($info_block, $SV_type, $alt, $start, $stop, \@format, \@tumour_parts, \@normal_parts, \@filter_reasons);
		}
		
		elsif ($type eq 'delly'){
			( $SV_length, $chr2, $t_SR, $t_PE, @filter_reasons ) = delly($info_block, $start, $stop, $SV_type, \@filter_reasons);
		}	
 	   	
		# chrom_filter($chr, $chr2, \@filter_reasons);
				
		$SVs{$id} = [ @fields[0..10], $SV_type, $SV_length, $stop, $chr2, $t_SR, $t_PE, \@filter_reasons ];
		
		$info{$id} = [ [@format], [@format_long], [@info_long], [@tumour_parts], [@normal_parts], [%sample_info], [%info_block] ];
			
	}	
	return (\%SVs, \%info);
}
		
sub lumpy {
	my ($info_block, $SV_type, $alt, $start, $stop, $format, $tumour_parts, $normal_parts, $filters) = @_;
	
	my @filter_reasons = @{ $filters };
	
	my ($SV_length) = $info_block =~ /SVLEN=(.*?);/;
		
			
		my ($t_PE, $t_SR, $c_PE, $c_SR);
		my ($t_DP, $c_DP);
		
		my $depth_flag = 0;

		for ( my $i = 0; $i <= $#{ $format }; $i++ ){
		
			if ( $format->[$i] eq 'PE' ) {
				$t_PE = $tumour_parts->[$i];
				$c_PE = $normal_parts->[$i];
			}
			
			elsif ( $format->[$i] eq 'SR' ){
				$t_SR = $tumour_parts->[$i];
				$c_SR = $normal_parts->[$i];
			}
				
			elsif ( $format->[$i] eq 'DP' ){
				$t_DP = $tumour_parts->[$i];
				$c_DP = $normal_parts->[$i];
				
				# Flag if either control or tumour has depth < 10 at site
				if ( $t_DP <= 10 or $c_DP <= 10 ){
					push @filter_reasons, 'tumour_depth=' . $t_DP;
					push @filter_reasons, 'control_depth=' . $c_DP;	
				}
			}		
		}
		
		my $tumour_read_support = ($t_PE + $t_SR);		
		my $control_read_support = ($c_PE + $c_SR);
		$control_read_support += 0.001;
		$tumour_read_support += 0.001;
		
		# Filters:
		if ($tumour_read_support <= 3) {
			push @filter_reasons, 'tumour_reads=' . $tumour_read_support;
		}
			
		# Filter if # tumour reads supporting var is less than 5 * control reads 
		# Or if there are more than 2 control reads 
		if ( ($control_read_support/$tumour_read_support) > 0.2 or $control_read_support > 2 ) {
			push @filter_reasons, 'control_reads=' . $control_read_support;
		}		
			
		my ($chr2) = 0;

		if ($SV_type eq 'BND'){
			$chr2 = $alt =~ s/[\[\]N]//g;
			($chr2, $stop) = $alt =~ /(.+)\:(\d+)/;
			$SV_length = $stop - $start;
		}
		
	return ($SV_length, $chr2, $t_SR, $t_PE, @filter_reasons);
}

sub delly {
	my ($info_block, $start, $stop, $SV_type, $filters) = @_;
	
	my @filter_reasons = @{ $filters };
	
		
	my ($SV_length) = ($stop - $start);
	
		my ($t_SR, $t_PE) = (0,0);
		
		if ($info_block =~ /;SR=(\d+);/){
	   		$t_SR = $1;
	   	}
	   
	   if ($info_block =~ /;PE=(\d+);/){
	   		$t_PE = $1;
	   }
 	  		
		if ($start > $stop){
			die "Start bigger than stop - shouldn't be here!!!!";
		}
		
		my ($chr2) = 0;
				
		if ($SV_type eq 'TRA'){
			($chr2) = $info_block =~ /CHR2=(.*?);/;
		}
	return ($SV_length, $chr2, $t_SR, $t_PE, @filter_reasons );
}

sub summarise_variants {
	my $SVs = shift;
	
	my ($dels, $dups, $trans, $invs, $filtered) = (0,0,0,0);
	
	my %support_by_chrom;
	
	my $read_support;
	
	my %filtered_sv;	
	
	for (keys %{ $SVs } ){
   	   
       my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $filters) = @{ $SVs->{$_} };
	   
	   my @filter_reasons = @{ $filters };
	   
		if (scalar @filter_reasons > 0){
			$filtered++;
		}
		
	   $read_support = ($SR + $PE);
	   $support_by_chrom{$id} = [$read_support, $sv_type, $chr ];
	  	   	   	   
	   $dels++ if $sv_type eq 'DEL';
	   $dups++ if $sv_type eq 'DUP';
	   $trans++ if $sv_type eq 'BND' or $sv_type eq 'TRA';
	   $invs++ if $sv_type eq 'INV';
	  
	}
	
	
	my ($max, @maxKeys);
	for (keys %support_by_chrom) {
	  $max = $support_by_chrom{$_}[0] unless $max;
	  if ($support_by_chrom{$_}[0] > $max) {
	      @maxKeys = ($_);
	      $max = $support_by_chrom{$_}[0];
	  } elsif ($support_by_chrom{$_}[0] == $max) {
	      push @maxKeys, $_;
	  }
	}

	print "\n";
	say "Variants where both break points are on chromosomes 2/3/4/X/Y:";
	say "$dels deletions";
	say "$dups duplications";
	say "$trans translocations";
	say "$invs inversions";
	
	print "\nSV with highest read support is:\nID: $_\nTYPE: $support_by_chrom{$_}[1]\nCHROM: $support_by_chrom{$_}[2]\nREADS: $support_by_chrom{$_}[0]\n" for @maxKeys;
	
}


sub get_variant {
	
	my ($id_lookup, $SVs, $info) = @_;
	
	if (not $info->{$id_lookup}){
		say "Couldn't find any variant with ID: '$id_lookup' in file. Abort";
		exit;
	}
	
	my (@format) 		= @{ $info->{$id_lookup}->[0]};
	my (@format_long) 	= @{ $info->{$id_lookup}->[1]};
	my (@info_long)		= @{ $info->{$id_lookup}->[2]};
	my (@tumour_parts) 	= @{ $info->{$id_lookup}->[3]};
	my (@normal_parts) 	= @{ $info->{$id_lookup}->[4]};
	my (%sample_info)	= @{ $info->{$id_lookup}->[5]};
	my (%info_block)	= @{ $info->{$id_lookup}->[6]};
	
			
	my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $filters) = @{ $SVs->{$id_lookup} };
    
	my @filter_reasons = @{ $filters };
	
	if (scalar @filter_reasons > 0){
	say "\n______________________________________________";	
	say "Variant '$id' filtered for the following reasons:";
	say "* $_" foreach @filter_reasons;
	say "______________________________________________\n";
	}
	
	say "ID:     $id";
	say "TYPE:   $sv_type";
	$chr2 ? say "CHROM1: $chr" : say "CHROM:  $chr";
	say "CHROM2: $chr2" if $chr2;	
	say "START:  $start";
	say "STOP:   $stop";
	say "LENGTH: $SV_length";
	say "PE:	$PE";
	say "SR:	$SR";
	say "QUAL:   $quality_score";
	say "FILT:   $filt";
	say "REF:    $ref";
	say "ALT:    $alt";
	
	say "____________________________________________________________________________________";
	printf " %-12s %-12s %-12s %-s\n", "INFO", "TUMOUR", "NORMAL", "EXPLAINER";
	say "____________________________________________________________________________________";
	 
	
	for (my $i = 0; $i <= $#format; $i++){
		printf " %-12s %-12s %-12s %-s \n", $format[$i], $tumour_parts[$i], $normal_parts[$i], $format_long[$i];
	}
	say "____________________________________________________________________________________";
	printf " %-12s %-12s %-s\n", "INFO", "VALUE", "EXPLAINER";
	say "____________________________________________________________________________________";
	
	for (sort keys %sample_info){
		printf " %-12s %-12s %-s \n", $_, $sample_info{$_}, $info_block{$_};
	}
	say "____________________________________________________________________________________";
	
	# say "INFO:   $info_block";
# 	say "FORMAT: $format_block";
# 	say "TUMOUR: $tumour_info_block";
# 	say "NORMAL: $normal_info_block";
		
}


sub chrom_filter {
	my $file = shift;
	
	my ($SVs, $info) = typer($file);
		
	my @keys = qw / 2L 2R 3L 3R 4 X Y /;
	my %chrom_filt;

	$chrom_filt{$_} = 1 for (@keys);
	
	for ( keys %{ $SVs } ){
    	my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $filter_reasons) = @{ $SVs->{$_} };
	
		if ($chr2 eq '0'){
			$chr2 =$chr;
		}
	
		if (not $chrom_filt{$chr} or not $chrom_filt{$chr2}){
			push @{ $filter_reasons }, 'chrom1=' . $chr . ';chrom2=' . $chr2;
			delete $SVs->{$id};
		}
	}
	return (\%$SVs, $info);
}



1;