#!/usr/bin/perl
#
use strict;
use warnings;

my $debug = '';

my $version = '
Version: 0.1.1 (2016-08-12)
Author: Adam Yongxin Ye @ CBI
';
my $usage = "Usage: $0 <.bam> <.bed>
Output:  STDOUT
    append 3 columns to <.bed>:
	read_count_for_each_signal(separated by \",\")
	forward_start(min,max)  reverse_start(min,max)
".$version;

if(@ARGV<2){
	die $usage;
}

my $bam_filename = shift(@ARGV);
my $bed_filename = shift(@ARGV);

my $bed_line;
my @F;
my ($bed_chr, $bed_start, $bed_end, $bed_strand);
my $strand_flag; 
my $command;
my $cmd_output;
my @sam_lines;
my @G;
my ($flag, $reverse, $duplicate);
my ($sam_chr, $sam_start);
my ($cigar, $cigar_array, $cigar_element, $cigar_end);
my (%main, $key, $value, @value_rearranged);
my ($read_name);
my (@starts);
open(BED, $bed_filename) or die "Error: cannot open the bed file $bed_filename\n";
while($bed_line = <BED>){
	chomp($bed_line);
	@F = split(/\t/, $bed_line);
	($bed_chr, $bed_start, $bed_end, $bed_strand) = @F[0,1,2,5]  ;   # bed format: 0-based, end-excluded 0-indexed  so 5 is strand SM -1/31/2024

	$bed_start++;   # convert to 1-based, end-included
	#print STDERR "[strand] $bed_strand \t" ;
	$strand_flag = "";
	#print STDERR ($bed_strand cmp '-'); 
	if(($bed_strand cmp '-') == 0){
		#print STDERR "[if was true ] \n";
		$strand_flag = "-f 16"; 
	}else{
		#print STDERR "[if was false ] \n "; 
		$strand_flag = "-F 16";
	}
	#print STDERR "[flag] $strand_flag\n"; 
	$command = "samtools view $strand_flag $bam_filename $bed_chr:$bed_start-$bed_end"; #add a -f command to get the strand that the reads aligned to 
	print STDERR "[CMD] $command\n";
	$cmd_output = `$command`;
	if($debug){
        	print STDERR "[DEBUG]\ttest debug\n";
        }
	
	@sam_lines = split(/[\n\r]+/, $cmd_output);
	%main = ();
	@starts = (
		[-1,-1],   # forward min,max_start
		[-1,-1],   # reverse min,max_start
	);
	foreach (@sam_lines){
		chomp;
		@G = split/\t/;
		$flag = $G[1];
		$reverse = $flag & 0x10;   # reverse
		$duplicate = $flag & 0x400;   # duplicate
		$reverse = 0 + ($reverse > 0);
		$duplicate = 0 + ($duplicate > 0);
		if($debug){
                        print STDERR "[DEBUG]\t$reverse\t$duplicate\n";
                }		
		$sam_chr = $G[2];
		$sam_start = $G[3];
		if($reverse){
			$cigar = $G[5];
#			print STDERR "[DEBUG] cigar = $cigar\n";
			$cigar_array = split_CIGAR($cigar);
#			print STDERR "[DEBUG] cigar_array = " . join(",", @$cigar_array) . "\n";
			foreach $cigar_element (@$cigar_array){
				$cigar_end = str_end($cigar_element, 1);
#				print STDERR "[DEBUG] cigar_element = $cigar_element\n";
#				print STDERR "[DEBUG] cigar_end = $cigar_end\n";
#				print STDERR "[DEBUG] str_trim($cigar_element, 1) = ".str_trim($cigar_element, 1)."\n";
				if($cigar_end eq "S" || $cigar_end eq "H" || $cigar_end eq "I" || $cigar_end eq "P"){
					# do nothing, but skip
				}elsif($cigar_end eq "M"|| $cigar_end eq "X" || $cigar_end eq "=" || $cigar_end eq "D"  || $cigar_end eq "N" ){
					$sam_start += str_trim($cigar_element, 1);
				}else{
					print STDERR "Warning: unrecognized cigar_end = $cigar_end\n";
				}
			}
			$sam_start--;   # end = start + len - 1
		}
		
		$key = $sam_chr . ":" . $sam_start . ":" . $reverse;
		if($debug){
			$read_name = $G[0];
			print STDERR "[DEBUG]\t$read_name\t$key\t$reverse\t$duplicate\n";
		}
		if(!exists($main{$key})){
			$main{$key} = [0, 0];   # non-duplicates, duplicates
		}
		if($duplicate){
			$main{$key}->[1]++;
		}else{
			$main{$key}->[0]++;
		}
		
		if($starts[$reverse]->[0] == -1 || $sam_start < $starts[$reverse]->[0]){
			$starts[$reverse]->[0] = $sam_start;
		}
		if($starts[$reverse]->[1] == -1 || $sam_start > $starts[$reverse]->[1]){
			$starts[$reverse]->[1] = $sam_start;
		}
	}
	
	@value_rearranged = ();
	foreach $key (sort {
				my $value_a = $main{$a};
				my $value_b = $main{$b};
				my $sum_a = $value_a->[0] + $value_a->[1];
				my $sum_b = $value_b->[0] + $value_b->[1];
				$sum_b <=> $sum_a || $value_b->[0] <=> $value_a->[0];
			} keys %main){
		$value = $main{$key};
#		push(@value_rearranged, ($value->[0] + $value->[1]) . "(" . $value->[0] .")");
		push(@value_rearranged, ($value->[0] + $value->[1]) );
	}
	if (@value_rearranged == 0){
		push(@value_rearranged, (0) );
	}
	print $bed_line . "\t" . join(",", @value_rearranged);
	print "\t" . join(",", @{$starts[0]}). "\t" . join(",", @{$starts[1]});
	print "\n";
}
close(BED);


sub split_CIGAR{
	my $cigar = shift(@_);
	my @ans = ();
	# CIGAR String \*|([0-9]+[MIDNSHPX=])+
	while($cigar =~ /^(\d+\w)(.*)$/){
		push(@ans, $1);
		$cigar = $2;
	}
	return [@ans];
}

sub str_end{
	my ($str, $len) = @_;
	my $strlen = length($str);
	return substr($str, $strlen-$len, $len);
}
sub str_trim{
	my ($str, $len) = @_;
	my $strlen = length($str);
	return substr($str, 0, $strlen-$len);
}

