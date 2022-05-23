#!/usr/bin/perl -w

use Getopt::Std;
 &getopts ("hl:L:F:N:t:T:");


if ($opt_h) {
print <<HELP;
STDIN and STDOUT

This script analyzes data in TAREAN ppm_Xmer_Y.csv file to:
	1) extract oligonucleotide sequences of specified length
	2) count for each oligonucleotide sequence a total coverage by k-mers that were used for consensus sequence reconstruction in the TAREAN pipeline
	3) calculate for each oligonucleotide sequence melting temperature (Tm = 81.5+16.6*(LOG10(0.39))+41*(GC/100)-500/length-0.62*50F)


Usage example:
1) sorted by 1. oligonucleotide length and 2. position (STD output of the script)
cat ppm_11mer_2.csv | PN_TAREAN_probe-Tm_calc_01.pl -l 50 -L 50 -F 50 -N 0.390 -t 52 -T 53

2) sorted by total coverage by k-mers
cat ppm_11mer_2.csv | PN_TAREAN_probe-Tm_calc_01.pl -l 50 -L 50 -F 50 -N 0.390 -t 52 -T 53 | sort -t\$'\\t' -k3,3n


Parameters (default values are in the square brackets):
-F [50] percentage of formamide 
-l [50] minimum oligonucleotide length, it should be 50-100 nt
-L [60] maximum oligonucleotide length, it should be 50-100 nt
-N [0.390] mM concentration of Na+ (0.390 mM = 2xSSC)
-t [0] minimum Tm 
-T [100] maximum Tm 

-h this help

HELP
exit;
}

if (defined $opt_l) {
	$min_kmer_length = $opt_l;
} else {
	$min_kmer_length = 50;
}


if (defined $opt_L) {
	$max_kmer_length = $opt_L;
} else {
	$max_kmer_length = 50;
}

if (defined $opt_F) {
	$formamide_perc = $opt_F;
} else {
	$formamide_perc = 50;
}

if (defined $opt_N) {
	$conc_Na = $opt_N;
} else {
	$conc_Na = 0.39;
}

if (defined $opt_t) {
	$min_Tm = $opt_t;
} else {
	$min_Tm = 0;
}

if (defined $opt_T) {
	$max_Tm = $opt_T;
} else {
	$max_Tm = 100;
}


$sequence_length = 0;

while ($line = <>) {
	chomp $line;
	if ($line =~ /^(\d+.*\d*)\t(\d+.*\d*)\t(\d+.*\d*)\t(\d+.*\d*)\t(\d+.*\d*)/) {
		# print "$line\n";
		@kmer_coverages = ($1,$2,$3,$4,$5);
		$max_cov = 0;
		for ($i=0; $i < @kmer_coverages; $i++) {
			if ($kmer_coverages[$i] > $max_cov) {
				$max_cov = $kmer_coverages[$i];
				if ($i == 0) {
					$nt = A;
				} elsif ($i == 1) {
					$nt = C;
				} elsif ($i == 2) {
					$nt = G;
				} elsif ($i == 3) {
					$nt = T;
				} elsif ($i == 4) {
					$nt = N;
				} else {
				 die "Wrong number of columns\n";
				}
			}
		}
		push (@max_coverages, $max_cov);
		push (@nts_with_max_coverage,$nt);
		$sequence_length ++;
	} elsif ($line =~ /^A\tC\tG\tT\t-/) {
		next;
	} elsif ($line =~ /^\s*$/) {
		next;
	} else {
		die "Unexpected structure of the line:\n$line\n";
	}
}



print "Position\tSequence\tCoverage_score\tLength\tGC_percentage\tTm\tNote\n";


$kmer_length = $min_kmer_length;

while ($kmer_length <= $max_kmer_length) {
	&extract_kmers;
	$kmer_length++;
}

sub extract_kmers {
	$position = 0;
	for ($i=0; $i<$sequence_length; $i++) {
		$position ++;
		$kmer_sequence = "";
		$kmer_total_coverage = 0;
		for ($k=0; $k < $kmer_length; $k++) {
			if ($i+$k < $sequence_length) {
				$index = $i+$k;
			} else {
				$index = $i+$k;
				while ($index >= $sequence_length) {
					$index = $index - $sequence_length;
				}
			}
			$kmer_sequence .=  $nts_with_max_coverage[$index];
			$kmer_total_coverage += $max_coverages[$index];
			#print "$index\n"
		}
	
		$GC_from_kmer_sequence = $kmer_sequence;
		$GC_from_kmer_sequence = $GC_from_kmer_sequence =~ s/[AT]//rg;
		$GC_count = length($GC_from_kmer_sequence);
		$GC_perc = $GC_count /$kmer_length * 100;
		$Tm = 81.5 + 16.6 * (log(0.39)/log(10)) + 41 * ($GC_perc / 100) - 500 / $kmer_length - 0.62 * $formamide_perc;
		if ($Tm >= $min_Tm and $Tm <= $max_Tm) {
			print "$position\t$kmer_sequence\t";
			printf "%0.5f", $kmer_total_coverage;
			print "\t$kmer_length\t";
			printf "%0.2f",$GC_perc;
			print "\t";
			printf "%0.2f",$Tm;
			if ($kmer_length < 50) {
				print "\tTm is unreliable, use other method of Tm calculation (e.g. dnaMATE)\n";
			} else {
				print "\t-\n";
			}
		}
	}
}





































