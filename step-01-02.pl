#!/usr/bin/perl

#####################################
#                                   #
#     Written by Kai Battenberg     #
#     Plant Sciences UC Davis       #
#                                   #
#####################################

use strict;
use warnings;

use Getopt::Long;
use Cwd;
use File::Path;

#####SCRIPT DESCRIPTION#####
#Script "step-01-02.pl" runs a quality check on all queries and generates a single-lined FASTA file with all queries.
#All output files will be in the working directory.
##########



#####Options#####
#setting the default values.
my $help = "";
my $query = ""; #file path to the query files.
my $q_seq_type = ""; #"DNA" or "AA".
my $q_clean = "YES"; #"YES" to clean the headers or "NO" to keep headers as is.
my $q_full_length = "NO"; #"YES" to only keep full-length sequences or "NO" to not have this screen.
my $q_qc_handle = "REPORT"; #"REMOVE" to automatically remove problematic sequences or "REPORT" to report and kill the process when such sequence is found.
my $spp_list = ""; #File with the species list.

#making the options into external arguments.
GetOptions (
	'help' => \$help,
	'query=s@' => \$query,
	'q_seq_type=s' => \$q_seq_type,
	'q_clean=s' => \$q_clean,
	'q_full_length=s' => \$q_full_length,
	'q_qc_handle=s' => \$q_qc_handle,
	'spp_list=s' => \$spp_list
);

#printing help
if ($help) {
	print "\n";
	print "\tstep-01-02.pl generates a single fasta file with all the queries formatted for the downstream processes of OrthoReD based on all the query sequence file(s) provided.\n\n";
	print "\t--help\tPrint this massage and die.\n\n";
	print "\t--query\tFile path to the directory with the query files.\n\n";
	print "\t--q_seq_type\tSequence type of the query. Values 'DNA' or 'AA' are accepted.\n\n";
	print "\t--q_clean\t'YES' to format the headers of the query sequences using the module <header_cleaner.pm> or 'NO' to leave the headers unchanged. (DEFAULT: YES)\n\n";
	print "\t--q_qc_handl\t'REPORT' to report and kill the process when a sequence with unclear reading frame is encountered or 'REMOVE' to automatically remove such sequence without killing the process. (DEFAULT: REPORT)\n\n";
	print "\t--q_full_length\t'YES' to only keep full-length sequences or 'NO' to not have this screen. (DEFAULT: NO)\n\n";
	print "\t--spp_list\tPath to the file with the lsit of species. See spp_list_EXAMPLE.txt for the proper format.\n";
	die "\n";
}

#checking for required options.
if (!$query) {
	die "USAGE: option --query is required.\n";
}
if ($q_seq_type !~ m/^DNA$/ && $q_seq_type !~ m/^AA$/) {
	die "USAGE: option --q_seq_type must be set as DNA or AA.\n";
}
if (!$spp_list) {
	die "USAGE: option --spp_list is required.\n";
}
##########



#####(Procedure-00)Setting up the environment.#####

###Identifying the working directory.
my $pwd = cwd();
my $pwd_ortho;
my @ls;
my @ls_ortho;

opendir my $dh, $pwd or die "Couldn't open dir '$pwd': $!";
@ls = readdir $dh;
closedir $dh;
###

###Checking all the required files###
#"OrthoReD_vXXXXXXXX" Directory with all scripts needed.
my $test_st = 1;
foreach my $file (@ls) {
	if ($file =~ m/OrthoReD/) {
		$test_st = 0;
		
		$pwd_ortho = $pwd."/".$file;
		opendir my $dh, $pwd_ortho or die "Couldn't open dir '$pwd_ortho': $!";
		@ls_ortho = readdir $dh;
		closedir $dh;
		
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-01-02.pl> requires folder <OrthoReD_vXXXXXXXX> to be in the working directory.\n";
}
else {
	$test_st = 1;
}

#"step-01-02.pl" script for this process.
foreach my $file (@ls_ortho) {
	if ($file =~ m/step-01-02.pl/) {
		$test_st = 0;
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-01-02.pl> requires file <step-01-02.pl> to be in folder <OrthoReD_vXXXXXXXX>.\n";
}
else {
	$test_st = 1;
}

#"OrthoReD_library.pm" Perl module with all the functions.
foreach my $file (@ls_ortho) {
	if ($file =~ m/OrthoReD_library.pm/) {
		$test_st = 0;
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-01-02.pl> requires file <OrthoReD_library.pm> to be in folder <OrthoReD_vXXXXXXXX>.\n";
}
else {
	$test_st = 1;
}

#"header_cleaner.pm" Perl module with header_cleaner function.
foreach my $file (@ls_ortho) {
	if ($file =~ m/header_cleaner.pm/) {
		$test_st = 0;
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-01-02.pl> requires file <OrthoReD_library.pm> to be in folder <OrthoReD_vXXXXXXXX>.\n";
}
else {
	$test_st = 1;
}
###

###setting up to use the libraries
use OrthoReD_vNFC::OrthoReD_library;
use OrthoReD_vNFC::header_cleaner;
###

##########



#####Programs required#####
#No external programs required.
##########



#####Generating required files#####
#Making the query directory
my $query_dir = "Step-01_Queries_RAW";
rmtree($query_dir);
mkdir $query_dir;

#Making the product directory
my $product_dir = "Step-02_QUERY";
rmtree($product_dir);
mkdir $product_dir;

#Making a log file
my $log = "step-01-02_LOG.txt";
open (LOG, ">$log") or die "cannot open $log.\n";
print "RUNNING SCRIPT: $pwd_ortho step-01-02.pl\n";
print LOG "RUNNING SCRIPT: $pwd_ortho step-01-02.pl\n";
print "(Procedure-00)Setting up the environment.\n";
print LOG "(Procedure-00)Setting up the environment.\n";

#Identify all the species included in the the species list.
my $full_name = 0;
my $short_name = 1;
my $acronym_name = 2;

my %spp1 = %{ &OrthoReD_library::species_lister ($spp_list, $acronym_name, $full_name) };
my %spp2 = %{ &OrthoReD_library::species_lister ($spp_list, $short_name, $full_name) };

print "Following species are included in the species list:\n";
print LOG "Following species were included in the species list:\n";
foreach my $sp (sort keys %spp1) {
	print "\t$spp1{$sp}\n";
	print LOG "\t$spp1{$sp}\n";
}

#Identifying all the raw query files and copy them into directory Step-01_Queries_RAW.
my @queries = @{ $query };
print "Following files are included in the query:\n";
print LOG "Following files were included in the query:\n";

foreach my $file (@queries) {
	print "\t$file\n";
	print LOG "\t$file\n";
	
	my @file_name = split (/\//, $file);
	my $new_file = $pwd."/".$query_dir."/".$file_name[-1];
	
	system "cp $file $new_file";
}

##########



#####(Procedure-01)Generate one quality-checked query file in single-line FASTA format.#####
print "(Procedure-01)Generate one quality-checked query file in single-line FASTA format.\n";
print LOG "(Procedure-01)Generate one quality-checked query file in single-line FASTA format.\n";

#Setting file names
my $procedure_01_dir = "01_QUERY_QCD";
my $qcd_file = "QUERY.fas";
my $removed_file = "REMOVED.fas";

if ($q_full_length =~ m/YES/) {
	print "Sequences that are not full-length CDS are removed.\n";
	print LOG "Sequences that are not full-length CDS were removed.\n";
}
elsif ($q_full_length =~ m/NO/) {
	print "Sequences with undetermined reading frame are removed.\n";
	print LOG "Sequences with undetermined reading frame were removed.\n";
}
else {
	die "ERROR: --q_full_length needs to be YES or NO.\n";
}

mkdir $procedure_01_dir;

$removed_file = $procedure_01_dir."/".$removed_file;
$qcd_file = $procedure_01_dir."/".$qcd_file;

#Generating one file with all queries.
my $query_files = join (" ", @queries);
system "cat $query_files > temp1.fas";

&OrthoReD_library::singler ("temp1.fas", "temp2.fas");

#Running quality control on sequences.
if ($q_clean =~ m/YES/) {
	open (TEMP, "<temp2.fas") or die "cannot open temp2.fas.\n";
	open (CLEAN, ">temp3.fas") or die "cannot open temp3.fas.\n";
	while (my $line = <TEMP>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		if ($line =~ m/^>/) {
			my $header = $line;
			my ($cleaned_header) = &header_cleaner::header_cleaner ($spp_list, $header);
			print CLEAN "$cleaned_header\n";
			next;
		}
		else {
			print CLEAN "$line\n";
			next;
		}
	}
	close (TEMP);
	close (CLEAN);
}
elsif ($q_clean =~ m/NO/) {
	system "cp temp2.fas temp3.fas";
}
else {
	die "ERROR: --q_clean needs to be YES or NO.\n"
}

my $value = "TRUE";
$value = &OrthoReD_library::qcer_file ("temp3.fas", $q_seq_type, $q_full_length);

if ($value =~ m/^TRUE$/) {
	system "cp temp3.fas $qcd_file";
	system "touch $removed_file";
	print "0\% of the sequences were removed.\n";
	print LOG "0\% of the sequences were removed.\n";
}
else {
	if ($q_qc_handle =~ m/^REPORT$/) {
		die "$value\n";
	}
	elsif ($q_qc_handle =~ m/^REMOVE$/) {
		my ($percent_removed) = &OrthoReD_library::separator ($spp_list, "temp3.fas", $qcd_file, $removed_file, $q_seq_type, $q_full_length);
		print "$percent_removed\% of the sequences were removed.\n";
		print LOG "$percent_removed\% of the sequences were removed.\n";
	}
	else {
		die "ERROR: --q_qc_handle needs to be REPORT, or REMOVE.\n";
	}
}

unlink ("temp1.fas");
unlink ("temp2.fas");
unlink ("temp3.fas");
##########



#####(Procedure-02)Formatting directories.#####
print "(Procedure-02)Formatting directories.\n";
print LOG "(Procedure-02)Formatting directories.\n";

#Formatting directories.
system "mv $procedure_01_dir $product_dir";
#rename ($procedure_01_dir, $product_dir);

print "SCRIPT COMPLETE\n";
print LOG "SCRIPT COMPLETE\n";
close (LOG);
system "mv $log $product_dir";
##########
__END__