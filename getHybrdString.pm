package getHybrdString;
#!/usr/bin/perl
use strict;
use warnings;


####################################################

sub get_hybrd_string {
	
	my $positive_fasta = shift or die;
	my $negative_fasta = shift or die;
	my $test_fasta = shift or die;
	my $positive_hybrd_string_fasta = shift or die;
	my $negative_hybrd_string_fasta = shift or die;
	my $test_hybrd_string_fasta = shift or die;

	
	#============ predict structure ===========
	
	my $positive_struct = "$positive_fasta" . "_struct";
	system("perl getStructure.pl $positive_fasta $positive_struct");

	my $negative_struct = "$negative_fasta" . "_struct";
	system("perl getStructure.pl $negative_fasta $negative_struct");
	
	my $test_struct = "$test_fasta" . "_struct";
	system("perl getStructure.pl $test_fasta $test_struct");
	
	
	#========== parse structure ==============
	
	my %hash_struct_positive = ();
	&parse_struct($positive_struct,\%hash_struct_positive);
	
	my %hash_struct_negative = ();
	&parse_struct($negative_struct,\%hash_struct_negative);
	
	my %hash_struct_test = ();
	&parse_struct($test_struct,\%hash_struct_test);
	
	
	#=========== construct hybrd strings =============
	
	
	open HYDSTPOS,">$positive_hybrd_string_fasta" or die;
	for my $id ( sort keys %hash_struct_positive ){
		
		my $seq = $hash_struct_positive{$id}{seq};
		my $struct = $hash_struct_positive{$id}{struct};
		
		my ($strA,$strB) = &reverse_align_bp($seq,$struct);
		
		my @tmp1 = split //,$strA;
		my @tmp2 = split //,$strB;
		
		$seq = "";
		for(my $i = 0; $i < length($strA); $i++){
			$seq .= "$tmp1[$i]$tmp2[$i]";
		}
		
		print HYDSTPOS ">$id\n$seq\n";
	}
	close HYDSTPOS;
	
	
	open HYDSTNEG,">$negative_hybrd_string_fasta" or die;
	for my $id ( sort keys %hash_struct_negative ){
		
		my $seq = $hash_struct_negative{$id}{seq};
		my $struct = $hash_struct_negative{$id}{struct};
		
		my ($strA,$strB) = &reverse_align_bp($seq,$struct);
		
		my @tmp1 = split //,$strA;
		my @tmp2 = split //,$strB;
		
		$seq = "";
		for(my $i = 0; $i < length($strA); $i++){
			$seq .= "$tmp1[$i]$tmp2[$i]";
		}
		
		print HYDSTNEG ">$id\n$seq\n";
	}
	close HYDSTNEG;
	
	
	open HYDSTEST,">$test_hybrd_string_fasta" or die;
	for my $id ( sort keys %hash_struct_test ){
		
		my $seq = $hash_struct_test{$id}{seq};
		my $struct = $hash_struct_test{$id}{struct};
		
		my ($strA,$strB) = &reverse_align_bp($seq,$struct);
		
		my @tmp1 = split //,$strA;
		my @tmp2 = split //,$strB;
		
		$seq = "";
		for(my $i = 0; $i < length($strA); $i++){
			$seq .= "$tmp1[$i]$tmp2[$i]";
		}
		
		print HYDSTEST ">$id\n$seq\n";
	}
	close HYDSTEST;
	
}

####################################################

sub parse_struct {
	
	my $infile_struct = shift or die;
	my $hash_struct = shift or die;
	
	unless(-e $infile_struct){
		print STDERR "Error:$infile_struct not exists\n";
		die;
	}
	
	open IN,"<$infile_struct" or die;
	
	my $id = "";
	my %hash_tmp = ();
	
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

sub is_numeric {

	use Scalar::Util qw(looks_like_number);

	my $v = shift;

	if( looks_like_number( $v ) ){  	
		return 1;
	} else {
		return 0;
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
	$lattice[0][0]{score}   = 0;
	$lattice[0][0]{pointer} = "none";
	
	for(my $j = 1; $j <= length($seqA); $j++) {
		$lattice[0][$j]{score}   = $GAP * $j;
		$lattice[0][$j]{pointer} = "left";
	}
	
	for (my $i = 1; $i <= length($seqB); $i++) {
		$lattice[$i][0]{score}   = $GAP * $i;
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



####################################################

sub reverse_align_bp {
	
	my $seq = shift or die;
	my $struct = shift or die;
	
	
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
	
	return($strA,$strB);
	
}

####################################################
1;


