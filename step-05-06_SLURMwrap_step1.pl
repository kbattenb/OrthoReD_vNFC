#!/usr/bin/perl

#####################################
#                                   #
#     Written by Kai Battenberg     #
#     Plant Sciences UC Davis       #
#                                   #
#####################################

#setting up to use the perl libraries
use strict;
use warnings;
use OrthoReD_vNFC::OrthoReD_library;
use OrthoReD_vNFC::header_cleaner;

#geting external variables
my $query = $ARGV[0];
my $database = $ARGV[1];

#Making directories
my $indir = "Step-05_INPUT";
my $inindir = $indir."/individualqueries";
my $outdir = "Step-06_ORTHOLOG";

system "mkdir -p $indir";
system "mkdir -p $inindir";
system "mkdir -p $outdir";

#Moving query
my @query = split (/\//, $query);
my $qname = $query[-1];
my $qfile = $indir."/".$qname;
system "cp -n $query $indir";
my @list = @{ &OrthoReD_library::fasta_splitter ($qfile, $inindir) };
my $list = join (",", @list);

#Moving database
my $dbname = "DATABASE";
my $dfolder = $indir."/".$dbname;
system "mkdir -p $dfolder";
system "cp -rn $database/* $dfolder";

print "$list";
