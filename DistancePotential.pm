package DistancePotential;
#!/usr/bin/perl
use strict;
use warnings;


#################
use getWordSet;
use getHybrdString;
#################


my $min_lng2;
my $max_lng2;
my $bin = 20;
my $rbin = 10;
my @array_kmers = ();
my @array_kmer_nmismatch =  ("3#0");


####################################################

sub extract_feature {
	
	my $fasta_positive = shift or die;
	my $fasta_negative = shift or die;
	my $fasta_test = shift or die;
	
	my $hash_feature_positive = shift or die;
	my $hash_feature_negative = shift or die;
	my $hash_feature_test = shift or die;
	
	
	
	# ========== Construct the set of kmers =============
	
	unless( @array_kmers ){
		print "\n======  Construct the set of kmers ======\n\n";
		my @array_character = qw(A C G U -);
		@array_kmers = getWordSet::get_mismatch_profile(\@array_kmer_nmismatch,\@array_character);
	}
	
	
	
	my $min_lng = 100;
	my $max_lng = 0;
	for my $kmer ( @array_kmers ){
		if( length($kmer) > $max_lng ){
			$max_lng = length($kmer);
		}
		if( length($kmer) < $min_lng ){
			$min_lng = length($kmer);
		}
	}
	$min_lng2 = 2*$min_lng;
	$max_lng2 = 2*$max_lng;
	
	print "\n====== Max & min length of kmers: $max_lng2\t$min_lng2 ========\n";
	
	
	# ========== Construct hybrid string =============
	
	
	print "\n======  Construct hybrid string ======\n\n";
	my $training_positive_fasta = "$fasta_positive.hybridString.fa";
	my $training_negative_fasta = "$fasta_negative.hybridString.fa";
	my $test_fasta = "$fasta_test.hybridString.fa";
	getHybrdString::get_hybrd_string($fasta_positive,$fasta_negative,$fasta_test,
	$training_positive_fasta,$training_negative_fasta,$test_fasta);
	
	
	# ========== Calculate statistical potentials =============

	print "\n======  Calculate statistical potentials ======\n\n";
	my %hash_potential = ();
	&kmer_pair_potentials($training_positive_fasta,$training_negative_fasta,\%hash_potential);
	
	my $potential_file = "$fasta_positive.potential";
	#&save_potential($potential_file,\%hash_potential);
	
	
	# ================ Extract features ================
	
	print "\n======  Extract statistical potential features ======\n\n";
	
	&get_feature($training_positive_fasta,\%hash_potential,$hash_feature_positive);
	
	&get_feature($training_negative_fasta,\%hash_potential,$hash_feature_negative);
	
	&get_feature($test_fasta,\%hash_potential,$hash_feature_test);
	
}




###############################################

sub get_feature{
	
	my $fasta = shift or die;
	my $hash_potential = shift or die;
	my $hash_feature = shift or die;
	
	my %hash_seq = ();
	&readfasta($fasta, \%hash_seq);
	
	
	for my $id ( sort keys %hash_seq ){
	
		my $seq = $hash_seq{$id};
		my %hash_seq_tmp = ();
		$hash_seq_tmp{$id} = $seq;
		
		# check sequence
		$seq =~ tr/[acgtTun\.]/[ACGUUUNN]/;	
		if( not ($seq =~/^([A|C|G|U|N|-]+)$/) ){
			next;
		}
		
		
		# get energe score
		my %hash_score = (); 
		&kmer_pair_frequency(\%hash_seq_tmp,
		'potential'=>$hash_potential,'score'=>\%hash_score);
		
		
		# output score feature
		for my $lng ( $min_lng2 .. $max_lng2 ){
			my $s = 0;
			for my $b ( 0 .. $bin ){
				if(exists $hash_score{$lng}{$b}){
					$s += $hash_score{$lng}{$b};			
				}
			}
			$hash_feature->{$id}{"Distance-lng$lng"} = $s;
		}
		
		#####################################
	}
}


####################################################

sub kmer_pair_frequency {
	
	my $hash_seq = shift or die;
	
	my $hash_dist;
	my $hash_potential;
	my $hash_score;
	my $hash_kmers;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/dist/i) {$hash_dist=shift @_}
		if ($argument=~/potential/i) {$hash_potential=shift @_}
		if ($argument=~/score/i) {$hash_score=shift @_}
		if ($argument=~/kmers/i) {$hash_kmers=shift @_}
	}
	
	
	
	# **************************************** #
	
	for my $id ( sort keys %{$hash_seq} ){
		
		my $iii = 0;
		my %lattice = ();
		my $seq = $hash_seq->{$id};
		for my $kmer ( @array_kmers ){
			my %coordinate = ();
			my $nmismatch = $kmer =~ tr/*/*/;  
			&matchPattern($nmismatch,$kmer,$seq,\%coordinate);
			for my $coord ( keys %coordinate ){
				$lattice{$iii}{"x"} = $coord;
				$lattice{$iii}{"base"} = $kmer;
				$iii++;
			}
		}
		
		# ==========================
		
		for(my $i = 0; $i < $iii; $i++) {
			for(my $j = $i + 1; $j < $iii; $j++) {
				
				################################################
				next if  length($lattice{$i}{"base"}) != length($lattice{$j}{"base"});
				################################################
				
				my $r = abs($lattice{$i}{"x"}-$lattice{$j}{"x"});				
				my $b = int( $r / $rbin );
				$b = $bin if $b > $bin;
							
				my $base_i = $lattice{$i}{"base"};
				my $base_j = $lattice{$j}{"base"};	
				
				($base_i, $base_j) = sort($base_i, $base_j);
				
				my $base_ij = "$base_i#$base_j";
				my $base_lng = length($base_ij) - 1;

				if( defined $hash_potential ){					
					if (exists $hash_potential->{$base_ij}{$b}){
						my $u = $hash_potential->{$base_ij}{$b};				
						if(defined $hash_score){
							unless(exists $hash_score->{$base_lng}{$b}){
								$hash_score->{$base_lng}{$b} = $u;
							}else{
								$hash_score->{$base_lng}{$b} += $u;
							}
						}
						if(defined $hash_kmers){
							unless(exists $hash_kmers->{$base_ij}{$b}){
								$hash_kmers->{$base_ij}{$b} = $u;
							}else{
								$hash_kmers->{$base_ij}{$b} += $u;
							}
						}
					} 	
				}elsif(defined $hash_dist){				
					unless(exists $hash_dist->{$base_ij}{$b}){
						$hash_dist->{$base_ij}{$b} = 1;
					} else {
						$hash_dist->{$base_ij}{$b}++;
					}	
				}
			}
		}
	}
}
	

sub kmer_pair_potentials {
	
	my $positive_fasta = shift or die;
	my $negative_fasta = shift or die;
	my $hash_potentials = shift or die;
	
	
	# ================= Positive ==================
	
	my %hash_seq = ();
	&readfasta($positive_fasta,\%hash_seq);
	
	my %hash_pos_dist = ();
	&kmer_pair_frequency(\%hash_seq,'dist'=>\%hash_pos_dist);
	
	
	# ================= Negative ==================
	
	%hash_seq = ();
	&readfasta($negative_fasta,\%hash_seq);
	
	my %hash_neg_dist = ();
	&kmer_pair_frequency(\%hash_seq,'dist'=>\%hash_neg_dist);
	
	
	# ==========================================
	

	for my $base_ij ( keys %hash_neg_dist ){
		for my $b ( 0 .. $bin ){	
			if ((exists $hash_pos_dist{$base_ij}{$b})  and (exists $hash_neg_dist{$base_ij}{$b})){
				my $r = -log ($hash_pos_dist{$base_ij}{$b} / $hash_neg_dist{$base_ij}{$b});
				unless(exists $hash_potentials->{$base_ij}{$b}){
					$hash_potentials->{$base_ij}{$b} = $r;
				} else {
					print STDERR "Error:($base_ij,$b) repeat\n";
					die;
				}	
			}
		}
	}	
}



sub save_potential {
	
	my $outfile = shift or die;
	my $hash_potentials = shift or die;
	
	open OUTPOT,">$outfile" or die;
	for my $base_ij ( sort keys %{$hash_potentials} ){
		for my $b ( 0 .. $bin ){
			if(exists $hash_potentials->{$base_ij}{$b}){
				my $r = $hash_potentials->{$base_ij}{$b};
				print OUTPOT "$base_ij\t$b\t$r\n";
			}
		}
	}
	close OUTPOT;	
	
}

####################################################

sub readfasta {
	
	my $usage = "< fasta file > < hash reference >";
	my $infile = shift or die $usage;
	my $hash_seq = shift or die $usage;
	
	unless(-e $infile){
		print STDERR "Error:$infile not exists";
		die;
	}
	
	open IN, $infile || die;
	
	my $c=0;
	my $seqId;
	while (defined (my $line=<IN>)) {
		chomp($line);
		next if $line=~/^\s*$/;
		if ($line =~/^>/) {
			$line =~s/\s*$//g;	
			$line =~s/^\s*//g;	
			$seqId = substr($line,1);
			$c++;
		} else {
			$line =~s/\s*//g;
			$hash_seq->{$seqId}.=$line;
		}
	}
	close IN;
	
	return $c;
	
}

####################################################

sub matchPattern {
	
	my $nmismatch = shift;
	my $pattern = shift or die;
	my $subjectSeq = shift or die;
	my $hash_pos = shift;
	
	
	my @pattern = split '', $pattern;
	my @subjectSeq = split '', $subjectSeq;

	my $lng = length( $pattern );
	my %err_count;
	for my $i ( 0 .. @subjectSeq - $lng ) {
		my $n_err = 0 ;
		for my $j ( 0 .. @pattern - 1 ) {
			$n_err++ if( $pattern[$j] ne $subjectSeq[$i+$j] );
		}
		next if ( $n_err > $nmismatch );
		$err_count{$i} = $n_err;
	}

	my $n = 0;
	foreach( keys %err_count ) {
		if( $err_count{$_} == $nmismatch ){
			my $value = substr( $subjectSeq, $_, $lng );
			if( defined $hash_pos ){
				$hash_pos->{ $_ } = $value;
			}	
			$n++;
		}
	}
	
	return $n;
}


####################################################
1;