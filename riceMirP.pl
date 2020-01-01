########################################
#!/usr/bin/perl
########################################
use strict;
use warnings;
use Cwd;
use Getopt::Long  qw(:config bundling);
use File::Basename;
########################################
use lib dirname(__FILE__);
########################################
use FeatureSelect;
use RandomForest;
use Performance;
use DistancePotential;
use PositionPotential;
use STRUCT;
########################################
my $is_feat_struct = 1;
my $is_feat_pair_potential = 1;
my $is_feat_position_potential = 1;
########################################
my $positive_fasta;
my $negative_fasta;
my $test_fasta;
my $resultTxt;
my $dir;
########################################

GetOptions( 	
	"P=s" => \$positive_fasta,
	"N=s" => \$negative_fasta,
	"T=s" => \$test_fasta,
	"O=s" => \$resultTxt,
	"D=s" => \$dir,
);

########################################

unless( defined $positive_fasta ){&usage(); exit;}
unless( defined $negative_fasta ){&usage(); exit;}
unless( defined $test_fasta ){&usage(); exit;}

my $current_dir = getcwd;
unless(defined $dir){$dir = "$current_dir/results";}
unless(-e $dir){mkdir($dir);}
unless(defined $resultTxt){$resultTxt = "$dir/results.txt";}

########################################

my $start_time = &start();

&train_predict($positive_fasta,$negative_fasta,$test_fasta,$resultTxt,$dir);

&end($start_time);


####################################################

sub train_predict {
	
	my $training_positive_fasta = shift or die;
	my $training_negative_fasta = shift or die;
	my $test_fasta = shift or die;
	my $outfile_out = shift or die;
	my $dir = shift or die;

	
	# ==========  Feature extraction ============

	my %hash_feature_positive = ();
	my %hash_feature_negative = ();
	my %hash_feature_test = ();
	
	
	if( $is_feat_struct ){
		STRUCT::extract_feature($training_positive_fasta,$training_negative_fasta,$test_fasta,
		\%hash_feature_positive,\%hash_feature_negative,\%hash_feature_test);
	}
	
	if( $is_feat_pair_potential ){
		DistancePotential::extract_feature($training_positive_fasta,$training_negative_fasta,$test_fasta,
		\%hash_feature_positive,\%hash_feature_negative,\%hash_feature_test);
	}
	
	if( $is_feat_position_potential ){
		PositionPotential::extract_feature($training_positive_fasta,$training_negative_fasta,$test_fasta,
		\%hash_feature_positive,\%hash_feature_negative,\%hash_feature_test);
	}
	
	
	# ==========  Feature selection ==========

	my $train_feature = "$dir/train_feature.txt";
	my $test_feature = "$dir/test_feature.txt";
	FeatureSelect::feature_selection(\%hash_feature_positive,\%hash_feature_negative,
	\%hash_feature_test,$train_feature,$test_feature,$dir);
	

	# ============= Traning ==============
	
	my $training_model = "$dir/model.RData";
	RandomForest::training($train_feature,$training_model,$dir);
	
	# ============= Testing ==============
	
	RandomForest::predict($test_feature,$training_model,$outfile_out,$dir);
	
}

####################################################

sub start {
	
	my $out = shift; 

	my($second, $minute, $hour, $dayOfMonth, $month, 
	$yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$second = "0$second" if($second =~ /^\d$/);
	my $sTime = "$hour:$minute:$second";
	my $stime = time;
	my $start_time = localtime;

	print "\n\nstarted: $start_time\n\n";
	if( defined $out ){
		print $out "started: $start_time\n";
	}
    
	return $stime;
}


####################################################

sub end {
	
	my $stime = shift or die;
	my $out = shift;

	my $etime = time - $stime;
	my ($second, $minute, $hour, $dayOfMonth, $month, 
	$yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$second = "0$second" if($second =~ /^\d$/);
	my $eTime = "$hour:$minute:$second";

	my $end_time = localtime;
    
	print "\n\nended: $end_time\n";
	print "total:", int($etime / 3600),"h:",
	int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";

	if( defined $out ){
		print $out "\n\nended: $end_time\n";
		print $out "total:", int($etime / 3600),"h:",
		int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";
	}

	print "\n\n============== End ============\n\n";
}


####################################################

sub usage {
	
my $usage = << "USAGE";

Program: $0
Contact: Yao Yuangen <yygen89\@163.com>

Usage:

	+++++++++++++++++++++++++++++++++++++++++
		
	Options:
	-N [file] : fasta file containing negative sequences
	-P [file] : fasta file containing positive sequences
	-T [file] : fasta file containing test sequences
	
	+++++++++++++++++++++++++++++++++++++++++
	
	$0 [options] -P <postive fasta> -N <negative fasta> -T <test fasta>

	+++++++++++++++++++++++++++++++++++++++++ 

USAGE
print $usage;

}

