package PositionPotential;
#!/usr/bin/perl
use strict;
use warnings;



my $n_segment = 20; 
my @bps = qw(-A -C -G -U AA AC AG AU CC CG CU GG GU UU);


#############################################################

sub extract_feature {
	
	my $training_positive_fasta = shift or die;
	my $training_negative_fasta = shift or die;
	my $test_fasta = shift or die;
	
	my $hash_feature_positive = shift or die;
	my $hash_feature_negative = shift or die;
	my $hash_feature_test = shift or die;
	
	
	
	# ========== Predict secondary structure =============
	
	print "\n========  Predict secondary structure =======\n\n";
	
	my $training_positive_struct = "$training_positive_fasta.struct";
	system("perl getStructure.pl $training_positive_fasta $training_positive_struct");

	my $training_negative_struct = "$training_negative_fasta.struct";
	system("perl getStructure.pl $training_negative_fasta $training_negative_struct");

	my $test_struct = "$test_fasta.struct";
	system("perl getStructure.pl $test_fasta $test_struct");
	

	# ========== Calculate statistical potentials =============
	
	print "\n========  Calculate statistical potentials =======\n\n";
	my %hash_potential = ();
	&construct_potential($training_positive_struct,$training_negative_struct,\%hash_potential);
	
	
	# ================ Extract features ================
	
	print "\n======  Extract statistical potential features ======\n\n";
	&get_feature($training_positive_struct,\%hash_potential,$hash_feature_positive);
	
	&get_feature($training_negative_struct,\%hash_potential,$hash_feature_negative);
	
	&get_feature($test_struct,\%hash_potential,$hash_feature_test);

}


sub get_feature {
	
	my $struct_file = shift or die;
	my $hash_potential = shift or die;
	my $hash_feature = shift or die;
	
	
	my %hash_struct = ();
	&parse_struct($struct_file,\%hash_struct);

	for my $id ( sort keys %hash_struct ){
	
		my $seq = $hash_struct{$id}{seq};
		my $struct = $hash_struct{$id}{struct};
	
		my %hash_rabp = (); 
		&reverse_align_bp($seq,$struct,\%hash_rabp);
		
		
		#################################################
		##########   Knowledge-based energy features    ##########
		#################################################

		for my $bin ( 0 .. $n_segment-1 ){
			my $s = 0;
			for my $bp ( @bps ){
				if( (exists $hash_potential->{$bp}{$bin}) && ($hash_rabp{$bp}{$bin} > 0) ){
					$s += $hash_potential->{$bp}{$bin};
				}
			}
			$hash_feature->{$id}{"Position-$bin"} = $s;
		}
		
		
		for my $bp ( @bps ){
			my $s = 0;
			for my $bin ( 0 .. $n_segment-1 ){
				if( (exists $hash_potential->{$bp}{$bin}) && ($hash_rabp{$bp}{$bin} > 0) ){
					$s += $hash_potential->{$bp}{$bin};
				}
			}
			$hash_feature->{$id}{"Position-$bp"} = $s;
		}
		
		#################################################
		
	}
	
}


sub construct_potential {
	
	my $training_positive_struct = shift or die;
	my $training_negative_struct = shift or die;
	my $hash_potentials = shift or die;
	
	
	# ================= Positive ==================
	
	my %hash_struct = ();
	&parse_struct($training_positive_struct,\%hash_struct);

	
	my %hash_pos_dist = ();
	for my $id ( sort keys %hash_struct ){
		my $seq = $hash_struct{$id}{seq};
		my $struct = $hash_struct{$id}{struct};
		&reverse_align_bp($seq,$struct,\%hash_pos_dist);	
	}
	
	
	# ================= Negative ==================
	
	%hash_struct = ();
	&parse_struct($training_negative_struct,\%hash_struct);

	my %hash_neg_dist = ();
	for my $id ( sort keys %hash_struct ){
		my $seq = $hash_struct{$id}{seq};
		my $struct = $hash_struct{$id}{struct};
		&reverse_align_bp($seq,$struct,\%hash_neg_dist);	
	}
	
	# ==========================================
	
	
	for my $bp ( @bps ){
		for my $bin ( 0 .. $n_segment-1 ){
			if ( $hash_pos_dist{$bp}{$bin} > 0  and $hash_neg_dist{$bp}{$bin} > 0 ){	
				my $r = -log($hash_pos_dist{$bp}{$bin} / $hash_neg_dist{$bp}{$bin});
				unless(exists $hash_potentials->{$bp}{$bin}){
					$hash_potentials->{$bp}{$bin} = &round($r,6);
				} else {
					print STDERR "Error:($bp,$bin) repeat\n";
					die;
				}
			}
		}
	}	
	
}




sub reverse_align_bp {
	
	my $seq = shift or die;
	my $struct = shift or die;
	my $hash_freq = shift or die;
	
	
	my $rev_seq = "";
	my $rev_struct = "";
	my @seqs = split ('',$seq);
	my @structs = split ('',$struct);

	# "reverse" sequence and struct.
	for( my $i = $#seqs; $i >= 0; $i-- ){
		$rev_seq .= $seqs[$i];	
		if ( $structs[$i] eq ")" ){
			$rev_struct .= "(";
		} else {
			$rev_struct .= $structs[$i];
		}
	}
	

	# change struct string into sequence string,
	# then excute Needleman Wunsch alignement.
	$struct =~ tr/\)/\(/;
	$struct =~ tr/\(\./AT/;
	$rev_struct =~ tr/\(\./AT/;
	($struct, $rev_struct) = &Needleman_Wunsch( $struct, $rev_struct );
	$struct =~ tr/AT/\(\./;
	$rev_struct =~ tr/AT/\(\./;


	# get the sequence corresponding to struct.
	@seqs = split('',$seq);
	@structs = split('',$struct);
	
	my $j = 0;
	my $result_seq = "";
	
	foreach my $i ( 0 .. $#structs ){
		if ( $structs[$i] eq "-" ){
			$result_seq .= "-";
		} else {
			$result_seq .= $seqs[$j];
			$j++;
		}
	}
	
	
	my @rev_seqs = split('',$rev_seq);
	my @rev_structs = split('',$rev_struct);
	
	$j = 0;
	my $result_rev_seq = "";
	
	foreach my $i ( 0 .. $#rev_structs ){
		if ( $rev_structs[$i] eq "-" ){
			$result_rev_seq .= "-";
		} else {
			$result_rev_seq .= $rev_seqs[$j];
			$j++;
		}
	}
	
	
	
	# ****************************************** #
	
	
	my $strA = uc ( $result_seq );
	my $strB = uc ( $result_rev_seq );
	$strA =~ tr/T/U/;
	$strB =~ tr/T/U/;


	if ( length ( $strA ) ne length( $strB ) ){
		print STDERR "Error:the length of two strings is different\n";
		die;
	}
	
	my $string_lng = length ( $strA );
	my $segment_lng = int ( $string_lng / $n_segment ) + 1;
	
	my @arrayA = split ('',$strA);
	my @arrayB = split ('',$strB);
	
	my %hash_rabp_tmp = ();
	my %hash_total = ();
	for(my $i = 0; $i <= $#arrayA; $i++ ){
		my $t = int ( $i / $segment_lng );

		my $chA = $arrayA[$i];
		my $chB = $arrayB[$i];
		my @tmp = ($chA,$chB);
		@tmp = sort (@tmp);
		my $bp = join('',@tmp);		
		unless ( exists $hash_rabp_tmp{$t}{$bp} ){
			$hash_rabp_tmp{$t}{$bp} = 1;
		} else {
			$hash_rabp_tmp{$t}{$bp}++;
		}
		unless ( exists $hash_total{$bp} ){
			$hash_total{$bp} = 1;
		} else {
			$hash_total{$bp}++;
		}
	}
	

	for my $i (0 ..$n_segment-1){
		for my $j ( @bps ){
			unless ( exists $hash_rabp_tmp{$i}{$j} ){
				$hash_rabp_tmp{$i}{$j} = 0;
			} 
			unless ( exists $hash_total{$j} ){
				$hash_total{$j} = 1;
			}
		}
	}
	
	
	# normalization 
	for my $i ( sort {$a<=>$b} keys %hash_rabp_tmp ){
		for my $j ( sort  keys %{ $hash_rabp_tmp{ $i } } ){
			$hash_freq->{$j}{$i} += $hash_rabp_tmp{ $i }{ $j } / $hash_total{ $j };
		}
	}
	
}


####################################################
####################################################
####################################################

sub parse_struct {
	
	my $infile_struct = shift or die;
	my $hash_struct	= shift or die;
	
	unless(-e $infile_struct){
		print STDERR "Error:$infile_struct not exists\n";
		die;
	}
	

	my $id = "";
	my %hash_tmp = ();
	open IN,"<$infile_struct" or die;
	while( defined(my $line = <IN> ) ){
		chomp( $line );
		next if $line =~ /^\s*$/;	
			
		if( $line =~/>/ ){
			$id = substr($line,1);
		} elsif ( $line =~ /[ACGUTN]+/i ) {
			
			$line = uc( $line );
			$hash_tmp{$id}{seq} = $line;
				
		} elsif ( $line =~ /\.|\(|\)/  ){				
			
			my @tmp = split /\s+\(/,$line;	
			$hash_tmp{$id}{struct} = $tmp[0];
			
			my $mfe = $tmp[1];
			$mfe =~ s/\(|\)//g;			
			$hash_tmp{$id}{mfe}= $mfe;
			
		}
	}
	close IN;
	
	
	for my $id ( keys %hash_tmp ){
		
		my $mfe = $hash_tmp{$id}{mfe};
		my $seq = $hash_tmp{$id}{seq};
		my $struct = $hash_tmp{$id}{struct};
		
		unless( $struct =~ /\(/ ) {
			next;
		}
		
		unless( $struct =~ /\)/ ) {
			next;
		}
		
		if (length($struct) != length($seq)){
			next;
		}
		
		unless( &is_numeric($mfe) ) {
			next;
		}
	
		$hash_struct->{$id}{mfe} = $mfe;
		$hash_struct->{$id}{seq} = $seq;
		$hash_struct->{$id}{struct} = $struct;
		
	}

}


####################################################

sub Needleman_Wunsch {
	
	my $seqA  = shift or die;
	my $seqB  = shift or die;

	
	my $MATCH =  1; 
	my $MISMATCH = -1;
	my $GAP = -1; 


	my @lattice = ();
	$lattice[0][0]{score} = 0;
	$lattice[0][0]{pointer} = "none";
	
	for(my $j = 1; $j <= length($seqA); $j++) {
		$lattice[0][$j]{score} = $GAP * $j;
		$lattice[0][$j]{pointer} = "left";
	}
	
	for (my $i = 1; $i <= length($seqB); $i++) {
		$lattice[$i][0]{score} = $GAP * $i;
		$lattice[$i][0]{pointer} = "up";
	}

	
	for(my $i = 1; $i <= length($seqB); $i++) {
		for(my $j = 1; $j <= length($seqA); $j++) {
			my ($diagonal_score, $left_score, $up_score);
   
			my $letterA = substr($seqA, $j-1, 1);
			my $letterB = substr($seqB, $i-1, 1);	 
								   
			if ($letterA eq $letterB) {
				$diagonal_score = $lattice[$i-1][$j-1]{score} + $MATCH;
			 } else {
				$diagonal_score = $lattice[$i-1][$j-1]{score} + $MISMATCH;
			 }

			
			$up_score   = $lattice[$i-1][$j]{score} + $GAP;
			$left_score = $lattice[$i][$j-1]{score} + $GAP;

			
			if ($diagonal_score >= $up_score) {
				if ($diagonal_score >= $left_score) {
					$lattice[$i][$j]{score}   = $diagonal_score;
					$lattice[$i][$j]{pointer} = "diagonal";
				 } else {
					$lattice[$i][$j]{score}   = $left_score;
					$lattice[$i][$j]{pointer} = "left";
				 }
			 } else {
				if ($up_score >= $left_score) {
					$lattice[$i][$j]{score}   = $up_score;
					$lattice[$i][$j]{pointer} = "up";
				 } else {
					$lattice[$i][$j]{score}   = $left_score;
					$lattice[$i][$j]{pointer} = "left";
				 }
			 }
		 }
	}


	my $alignA = "";
	my $alignB = "";

	my $j = length($seqA);
	my $i = length($seqB);

	while ( 1 ) {
		
		last if $lattice[$i][$j]{pointer} eq "none"; 

		if ($lattice[$i][$j]{pointer} eq "diagonal") {
			$alignA .= substr($seqA, $j-1, 1);
			$alignB .= substr($seqB, $i-1, 1);
			$i--;
			$j--;
		 } elsif ($lattice[$i][$j]{pointer} eq "left") {
			$alignA .= substr($seqA, $j-1, 1);
			$alignB .= "-";
			$j--;
		 } elsif ($lattice[$i][$j]{pointer} eq "up") {
			$alignA .= "-";
			$alignB .= substr($seqB, $i-1, 1);
			$i--;
		 }	
		 
	}

	$alignA = reverse $alignA;
	$alignB = reverse $alignB;
	
	return ($alignA,$alignB);
}


sub is_numeric {
	
	use Scalar::Util qw(looks_like_number);
	
	my $v = shift;
	
	if( looks_like_number( $v ) ){  	
		return 1;
	} else {
		return 0;
	}
}

sub round {
	
	my $usage = "<  numerical value > < the numbers after the decimal point >";
	my $val = shift;
	my $col = shift;

	unless(defined $val){
		die $usage;
	}

	unless(defined $col){
		die $usage;
	}

	unless( &is_numeric($val) ){
		print STDERR "Error:$val not a numeric";
		die;
	}

	unless( &is_numeric($col) ){
		print STDERR "Error:$col not a numeric";
		die;
	}

	my $r = 10 ** $col;
	my $a = ($val > 0) ? 0.5 : -0.5;

	return int($val * $r + $a) / $r;
}

####################################################
1;