#!/usr/bin/perl -w
use Getopt::Std;
use File::Glob;

&getopts("hd:");

if ($opt_h) {
print <<HELP;
This script extracts and codes reads from TAREAN cluster directories (dir_CL0001, dir_CL0002 ...) in the TAREAN output

The structure of TAREAN cluster directories is as follows:
.../seqclust/clustering/clusters/dir_CL0001
.../seqclust/clustering/clusters/dir_CL0002
.../seqclust/clustering/clusters/dir_CL0003


Mandatory parameters:
-d path to the TAREAN "clusters" directory (.../seqclust/clustering/clusters)

Optional parameters:
-h this help

Example of usage:
PN_extract_and_code_reads_from_TAREAN_CLs_02.pl -d /mnt/ceph/454_data/workshop/2022/pavel/Protocol_3-TAREAN/seqclust/clustering/clusters > output.fasta

HELP
exit;
}

if (defined $opt_d) {
	$clusters_directory = $opt_d;
} else {
	die "-d was not specified";
}



@list_of_directories = grep {-d} glob ("$clusters_directory/*");

foreach $directory (@list_of_directories) {
	open (READS, "<$directory/reads.fasta");
	if ($directory =~ /^\S+\/(dir_CL\d+)$/) {
		$CL_directory = $1;
	} else {
		die "Unexpected structure of TAREAN directories:\n$directory\n"
	}
	while ($line = <READS>) {
		if ($line =~ /^>(\S+)/) {
			print ">$CL_directory"."_"."$1\n";
		} else { 
			print "$line";
		}
	}
	close READS;
	
}

