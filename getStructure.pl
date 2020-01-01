#!/usr/bin/perl
use strict;
use warnings;

###########################################
use RNA; # Vienna RNA package perl interface
###########################################


my $infile_fasta = shift or die;
my $outfile_struct = shift or die;

my %hash_seq = ();
&readfasta($infile_fasta,\%hash_seq);

open OUT,">$outfile_struct" or die;
for my $id ( keys %hash_seq ){
	my $seq = $hash_seq{$id};
	$seq = uc ( $seq );
	$seq =~ s/T/U/ig;
	$seq =~ s/^\s*//g;
	$seq =~ s/\s*$//g;
	unless( not ($seq =~/^([A|C|G|U|N]+)$/) ){
		my ($struct,$mfe) = RNA::fold($seq);
		$mfe = &round($mfe,4);
		print OUT ">$id\n$seq\n$struct ($mfe)\n";
	}	
}
close OUT;


###########################################

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

###########################################

sub round {
	
	my $val = shift;
	my $col = shift;

	my $r = 10 ** $col;
	my $a = ($val > 0) ? 0.5 : -0.5;

	return int($val * $r + $a) / $r;
}

###########################################