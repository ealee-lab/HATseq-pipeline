#!/usr/bin/perl

use strict;
use warnings;

my (@F, $i);
my $last_peak_name = "";
my $peak_name;
my %peak_position = ();
my $read_start_position;
my %peak_read_start_position = ();
my %peak_read_count = ();
while(<STDIN>){
	chomp;
	@F = split/\t/;
	
#	for($i=0; $i<@F; $i++){
#		print join("\t", $i, $F[$i])."\n";
#	}
	
	$peak_name = $F[9];
	$peak_position{$peak_name} = [@F[6..11]];
	
	if($last_peak_name ne "" && $peak_name ne $last_peak_name){
		print_peak($last_peak_name);
		$last_peak_name = "";
	}
	$last_peak_name = $peak_name;
	
	$read_start_position = $F[1] + 1;
	if($F[5] eq "-"){
		$read_start_position = $F[2];
	}
	
	if(!exists($peak_read_start_position{$peak_name})){
		$peak_read_start_position{$peak_name} = {};
	}
	$peak_read_start_position{$peak_name}->{$read_start_position} ++;
	$peak_read_count{$peak_name} ++;
}

if($last_peak_name ne ""){
	print_peak($last_peak_name);
	$last_peak_name = "";
}

#foreach $peak_name (sort {
#	$peak_position{$a}->[0] <=> $peak_position{$b}->[0] || $peak_position{$a}->[1] <=> $peak_position{$b}->[1] || $peak_position{$a}->[2] <=> $peak_position{$b}->[2];
#} keys %peak_read_start_position){
#	print join("\t", @{$peak_position{$peak_name}}, $peak_read_count{$peak_name}, scalar(keys %{$peak_read_start_position{$peak_name}}))."\n";
#}

0;


sub print_peak{
	my ($peak_name) = @_;
	print join("\t", @{$peak_position{$peak_name}}, $peak_read_count{$peak_name}, scalar(keys %{$peak_read_start_position{$peak_name}}))."\n";
}

