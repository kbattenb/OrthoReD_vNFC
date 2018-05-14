###############################################################################
#                                                                             #
#                           OrthoReD Perl module                              #
#                                                                             #
#                        Written by Kai Battenberg                            #
#                         Plant Sciences UC Davis                             #
#                                                                             #
#                     Copyright (C) 2017 Kai Battenberg                       #
#                                                                             #
#                                                                             #
#    Orthologfinder is free software. You can redistribute it and/or modify   #
#    it under the terms of the GNU General Public License (version 3 or any   #
#    later version) as published by the Free Software Foundation. The copy    #
#    of the GNU General Public License (version 3) is attached to the copy    #
#    of Orthologfinder, but if it somehow was not, see the following link     #
#    <http://www.gnu.org/licenses/gpl.txt> for more details.                  #
#                                                                             #
#    Orthologfinder is distributed as is without any warrenty. See the GNU    #
#    General Public License for more details.                                 #
#                                                                             #
###############################################################################

use strict;
use warnings;

use Getopt::Long;
use Cwd;
use File::Path;

package OrthoReD_library;

#####MODULES#####

#Given a single-lined AA FASTA file, tests if all 20 amino acids are present or not.
sub aminotester {
	if (@_ != 1) {
		die "USAGE: module aminotester requires one argument: <\$in_file>.\n";
	}
	
	#Taking in the arguments
	my ($in_file) = @_;
	
	#test results
	my $result = "PASS";
	
	#generating one string with all characters that are used
	my $tseq = `sed -n 'n;p' $in_file`;
	$tseq =~ s/\n//sig;
	$tseq =~ s/\r//sig;
	$tseq =~ s/\-//g;
	my @tseq = split (//, $tseq);
	
	my %seen = ();
	my @uniq;
	foreach my $item (@tseq) {
		push(@uniq, $item) unless $seen{$item}++;
	}
	@tseq = sort (@uniq);
	$tseq = join ("", @tseq);
	
	#testing if all 20 amino acids are present
	my $ref = "ACDEFGHIKLMNPQRSTVWY";
	my @ref = split (//, $ref);
	foreach my $aa (sort @ref) {
		if ($tseq =~ m/$aa/) {
			next;
		}
		else {
			$result = "FAIL";
			last;
		}
	}
	
	return ($result);
}

#Given a single-lined AA FASTA file and a database, runs AB-BLASTP for each sequence against the database.
sub blaster_AB {
	if (@_ != 4) {
		die "USAGE: module blaster_AB requires four arguments: <\$query> <\$db> <\$w> <\$id>.\n";
	}
	
	#Checking if BLAST is running.
	my $blast_check = `which xdformat`;
	if (!$blast_check) {
		die "USAGE: module blaster_AB requires AB-BLAST to be running.\n";
	}
	
	#Taking in the arguments
	my ($query, $db, $w, $id) = @_;
	
	#BLASTP conditions
	my $e = "1e-3";
	my $t = 1000;
	my $bv = 100000;
	
	#out_file
	my $out = "BLASTER_01_".$id.".txt";
	#temp_file
	my $temp = "BLASTER_02_".$id.".txt";
	
	#Begin BLAST.
	system "blastp $db $query mformat=2 wordmask=seg matrix=BLOSUM62 E=$e W=$w T=$t B=$bv V=$bv postsw warnings O=$temp >/dev/null";
	system "cut -f 1,2,3,6,7,11 $temp > $out";
	system "rm $temp";
	
	return ($out);
}

#Given a single-lined AA FASTA file and a database, runs NCBI-BLASTP for each sequence against the database.
sub blaster_NCBI {
	if (@_ != 4) {
		die "USAGE: module blaster_NCBI requires four arguments: <\$query> <\$masked_db> <\$id> <\$threads>.\n";
	}
	
	#Checking if BLAST is running.
	my $blast_check = `which segmasker`;
	if (!$blast_check) {
		die "USAGE: module blaster_NCBI requires NCBI-BLAST to be running.\n";
	}
	
	#Taking in the arguments
	my ($query, $masked_db, $id, $threads) = @_;
	
	#out_file
	my $out = "BLASTER_01_".$id.".txt";
	
	#BLAST conditions
	my $e_val = "1e-03";
	my $max_seqs = "25000";
	
	#BLAST command.
	my $blast_type = "blastp";
	my $outfmt = "6 qseqid sseqid evalue pident qcovs length";
	
	#mask id
	my $maskid = 21;
	my $infofile = $masked_db."_info.txt";
	my $subject;
	my $maskcheck = `[ -f $infofile ] && echo 'present' || echo 'absent'`;
	if ($maskcheck =~ m/present/) {
		$maskcheck = `grep -A1 '^Algorithm' $infofile | tail -n 1 | cut -d' ' -f1`;
		$maskcheck =~ s/\r//sig;
		$maskcheck =~ s/\n//sig;
	}
	else {
		$subject = `ls $masked_db*fas`;
		$subject =~ s/\r//sig;
		$subject =~ s/\n//sig;
	}
	
	#temp file
	my $temp = "BLASTER_01_temp_".$id."\.fas";
	
	system "touch $out";
	open (QUERY, "<$query") or die "cannot open $query.\n";
	while ($_ = <QUERY>) {
		if ($_ =~ m/^>/) {
			#Header
			my $header = $_;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			
			#Sequence
			my $seq = <QUERY>;
			$seq =~ s/\r//sig;
			$seq =~ s/\n//sig;
			
			open (TEMP_AA, ">$temp") or die "cannnot open $temp\n";
			print TEMP_AA "$header\n";
			print TEMP_AA "$seq\n";
			close (TEMP_AA);
			
			if ($maskcheck =~ m/^$maskid$/ && $maskid =~ m/^$maskcheck$/) {
				system "$blast_type -query $temp -db $masked_db -db_soft_mask $maskid -use_sw_tback -evalue $e_val -max_target_seqs $max_seqs -outfmt '$outfmt' -num_threads $threads >>$out";
			}
			elsif ($maskcheck !~ m/absent/) {
				system "$blast_type -query $temp -db $masked_db -use_sw_tback -evalue $e_val -max_target_seqs $max_seqs -outfmt '$outfmt' -num_threads $threads >>$out";
			}
			else {
				system "$blast_type -query $temp -subject $subject -use_sw_tback -evalue $e_val -max_target_seqs $max_seqs -outfmt '$outfmt' -num_threads $threads >>$out";
			}
			unlink ($temp);
		}
		else {
			next;
		}
	}
	close (QUERY);
	
	return ($out);
}

#Given a single-lined AA FASTA file and a database, runs SWIPE for each sequence against the database.
sub blaster_SWIPE {
	if (@_ != 4) {
		die "USAGE: module blaster_SWIPE requires four arguments: <\$query> <\$masked_db> <\$id> <\$threads>.\n";
	}
	
	#Checking if SWIPE is running.
	my $swipe_check = `which swipe`;
	if (!$swipe_check) {
		die "USAGE: module blaster_SWIPE requires SWIPE to be running.\n";
	}
	my $blast_check = `which blastdbcmd`;
	if (!$blast_check) {
		die "USAGE: module blaster_SWIPE requires NCBI-BLAST to be running.\n";
	}
	
	#Taking in the arguments
	my ($query, $masked_db, $id, $threads) = @_;
	
	#out_file
	my $out = "BLASTER_01_".$id.".txt";
	
	#SWIPE conditions
	my $e_val = "1e-03";
	
	#SWIPE command.
	my $blast_type = "blastp";
	my $outfmt = "8";
	
	#temp file
	my $temp1 = "BLASTER_01_temp1_".$id."\.fas";
	my $temp2 = "BLASTER_01_temp2_".$id."\.fas";
	
	system "touch $temp2";
	open (QUERY, "<$query") or die "cannot open $query.\n";
	while ($_ = <QUERY>) {
		if ($_ =~ m/^>/) {
			#Header
			my $header = $_;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			
			#Sequence
			my $seq = <QUERY>;
			$seq =~ s/\r//sig;
			$seq =~ s/\n//sig;
			
			open (TEMP_AA, ">$temp1") or die "cannnot open $temp1\n";
			print TEMP_AA "$header\n";
			print TEMP_AA "$seq\n";
			close (TEMP_AA);
			
			system "swipe --db=$masked_db --query=$temp1 --evalue=$e_val --num_threads=$threads --outfmt=$outfmt --symtype=$blast_type >>$temp2";
			unlink ($temp1);
		}
		else {
			next;
		}
	}
	close (QUERY);
	
	#get conversion list
	my $conv = `blastdbcmd -db $masked_db -entry all -outfmt "%o %t" 2>/dev/null`;
	$conv =~ s/\r$//;
	$conv =~ s/\n$//;
	my @conv = split (/\n/, $conv);
	my %conv;
	foreach my $element (@conv) {
		my @element = split (/ /, $element);
		$element[0] = "gnl\|BL_ORD_ID\|".$element[0];
		$conv{$element[0]} = $element[1];
	}
	
	#format temp2 file
	open (TEMP2, "<$temp2") or die "cannot open $temp2.\n";
	open (OUT, ">$out") or die "cannot open $out.\n";
	while (my $line = <TEMP2>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		my @line = split (/\t/, $line);
		
		print OUT "$line[0]\t$conv{$line[1]}\t$line[10]\t$line[11]\t$line[3]\t$line[2]\n";
	}
	close (TEMP2);
	close (OUT);
	unlink ($temp2);
	
	return ($out);
}

#Reads a tabulated BLAST output file and parses out hits hits according to the parameters that are set.
sub blast_parser {
	if (@_ != 6) {
		die "USAGE: module blast_parser requires six arguments: <\$in_file> <\$out_file> <\$eval_threshold> <\$identity_threshold> <\$length_threshold> <\$sander_schneider>.\n";
	}
	
	my ($in_file, $out_file, $eval_threshold, $identity_threshold, $length_threshold, $sander_schneider) = @_;
	
	my @eval_threshold = split (/e-/, $eval_threshold);
	$eval_threshold = 1 / (10 ** $eval_threshold[1]);
	
	open (IN_FILE, "<$in_file") or die "cannot open $in_file.\n";
	open (OUT_FILE, ">$out_file") or die "cannot open $out_file.\n";
	while (my $line = <IN_FILE>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		
		my @line = split (/\t/, $line);
		
		#setting the velues
		my $eval = 1;
		my $identity = 0;
		my $coverage = 0;
		my $length = 0;
		my $sander_schneider_test = "FAIL";
		
		#reading query and subject
		my $query = $line[0];
		my $subject = $line[1];
		
		#reading eval
		$eval = $line[2];
		if ($eval =~ m/e/) {
			my @value = split (/e-/,$eval);
			$eval = $value[0] / (10 ** $value[1]);
		}
		elsif ($eval =~ m/^0.$/) {
			$eval = 0;
		}
		
		#reading alignment length
		$length = $line[3];
		
		#reading %identity
		$identity = $line[4];
		
		#implementing Sander and Schneider (1991).
		if ($sander_schneider =~ m/YES/) {
			my $sander_schneider_threshold = 0;
			if ($length < 10) {
				$sander_schneider_threshold = 101;
			}
			elsif ($length >= 10 and $length <= 80) {
				$sander_schneider_threshold = 290.15 * ($length ** -0.562)
			}
			else {
				$sander_schneider_threshold = 25;
			}
			if ($identity >= $sander_schneider_threshold) {
				$sander_schneider_test = "PASS";
			}
		}
		
		#printing out lines
		if ($eval <= $eval_threshold && $identity > $identity_threshold && $length > $length_threshold) {
			if ($sander_schneider =~ m/YES/ && $sander_schneider_test =~ m/FAIL/) {
				next;
			}
			else {
				print OUT_FILE "$query\t$subject\t$eval\t$length\t$identity\n";
			}
		}
	}
	close (IN_FILE);
	close (OUT_FILE);
}

#Reads a tabulated BLAST output file and sorts the hits and query.
sub blast_organizer {
	if (@_ != 2) {
		die "USAGE: module blast_organizer requires two arguments: <\$in_file> <\$out_file>.\n";
	}
	
	my ($in_file, $out_file) = @_;
	
	open (IN, "<$in_file") or die "cannot open $in_file.\n";
	open (OUT, ">$out_file") or die "cannot open $out_file.\n";
	while (my $hit = <IN>) {
		$hit =~ s/\n//sig;
		$hit =~ s/\r//sig;
		
		my @hit = split (/\t/, $hit);
		my @qs = ($hit[0], $hit[1]);
		@qs = sort @qs;
		$hit[0] = $qs[0];
		$hit[1] = $qs[1];
		
		$hit = join ("\t", @hit);
		
		print OUT "$hit\n";
	}
	close (IN);
	close (OUT);
}

#Given a tree with branch lengths, removes the branch lengths
sub branch_length_stripper {
	if (@_ != 1) {
		die "USAGE: module branch_length_stripper requires one arguments: <\$tree>.\n";
	}
	
	my ($tree) = @_;
	
	my $out_tree = `nw_topology $tree`;
	$out_tree =~ s/\r//sig;
	$out_tree =~ s/\n//sig;
	$out_tree =~ s/\;$//;
	return ($out_tree);
	
}

#Reads single-lined DAN/AA FASTA file and generates a masking data and masked database in both DNA and AA.
sub database_maker_AB {
	if (@_ != 3 && @_ != 4) {
		die "USAGE: module database_maker_AB requires three or four arguments: <\$database> <\$sequence_type> <\$database_name> <\$path>.\n";
	}
	
	my $path = "";
	if (@_ == 3) {
		my ($database, $sequence_type, $database_name) = @_;
	}
	if (@_ == 4) {
		my $database = $_[0];
		my $sequence_type = $_[1];
		my $database_name = $_[2];
		$path = $_[3];
	}
	
	#
	#Checking if BLAST is running.
	#
	my $command = $path."xdformat";
	my $blast_check = `which $command`;
	if (!$blast_check) {
		die "USAGE: module database_maker_AB requires AB-BLAST to be running.\n";
	}
	
	#
	#Taking in the arguments.
	#
	my ($database, $sequence_type, $database_name) = @_;
	
	if ($sequence_type =~ m/^DNA$/) {
		system "$command -n -I -o $database_name $database";
	}
	if ($sequence_type =~ m/^AA$/) {
		system "$command -p -I -o $database_name $database";
	}
	mkdir ($database_name);
	system "mv $database_name\.* $database_name";
}

#Reads single-lined DAN/AA FASTA file and generates a masking data and masked database in both DNA and AA.
sub database_maker_NCBI {
	if (@_ != 3 && @_ != 4) {
		die "USAGE: module database_maker_NCBI requires three or four arguments: <\$database> <\$sequence_type> <\$database_name> <\$path>.\n";
	}
	
	my $database;
	my $sequence_type;
	my $database_name;
	my $path = "";
	if (@_ == 3) {
		$database = $_[0];
		$sequence_type = $_[1];
		$database_name = $_[2];
	}
	if (@_ == 4) {
		$database = $_[0];
		$sequence_type = $_[1];
		$database_name = $_[2];
		$path = $_[3];
	}
	
	#
	#Checking if BLAST is running.
	#
	my $test = $path."segmasker";
	my $blast_check = `which $test`;
	if (!$blast_check) {
		die "USAGE: module database_maker requires NCBI-BLAST to be running.\n";
	}
	
	my $masking_file = $database_name."_MASKING";
	my $database_log = $database_name."_log.txt";
	my $database_info = $database_name."_info.txt";
	
	my $cmd1;
	my $seqtype;
	if ($sequence_type =~ m/^DNA$/) {
		$cmd1 = $path."dustmasker";
		$seqtype = "nucl";
	}
	elsif ($sequence_type =~ m/^AA$/) {
		$cmd1 = $path."segmasker";
		$seqtype = "prot";
	}
	my $cmd2 = $path."makeblastdb";
	my $cmd3 = $path."blastdbcmd";
	
	system "$cmd1 -in $database -infmt fasta -outfmt maskinfo_asn1_bin -out $masking_file\.asnb";
	system "$cmd2 -in $database -dbtype $seqtype -input_type fasta -mask_data $masking_file\.asnb -title $database_name -out $database_name -logfile $database_log";
	system "$cmd3 -db $database_name -info -out $database_info";
	
	mkdir ($database_name);
	system "mv $database_log $database_name";
	system "mv $database_info $database_name";
	system "mv $masking_file\.asnb $database_name";
	system "mv $database_name\.* $database_name";
}

#Reads a tree file and generates the distance matrix of all tips.
sub distance_finder {
	if (@_ != 2) {
		die "USAGE: module distance_finder requires one arguments: <\$tre_file> <\$out_file>.\n";
	}
	
	my ($tre_file, $out_file) = @_;
	
	my $tree = `head -n 1 $tre_file`;
	$tree =~ s/\r//sig;
	$tree =~ s/\n//sig;
	
	#record the distance of parentheses
	my %dis;
	my $par = -1;
	my @current;
	while($tree =~ /./g) {
		if ($& eq '(') {
			$par ++;
			next if $par == 0;
			$current[$#current+1] = $par;
		}
		elsif ($& eq ')') {
			(my $tem) = $' =~ /:(\d+\.\d+|\d+)/;
			next if $#current == -1;
			$dis{'node_'.$current[$#current]} = $tem;
			pop @current;
		}
	}
	
	#record the distance of leaves
	my @order;
	while ($tree =~ /([^\(\):,]+):(\d+\.\d+|\d+)/g) {
		$dis{$1} = $2;
		$order[$#order+1] = $1;
	}
	
	#record parents of leaves
	my %pare;
	@current = ();
	$par = -1;
	while ($tree =~ /(\(|\)|([^\(\):,]+):)/g) {
		if ($& eq '(') {
			$par ++;
			if ($par == 0) {
				next;
			}
			$current[$#current+1] = $par;
		}
		elsif ($& eq ')') {
			pop @current;
		}
		else {
			map {$pare{$2}{$_} = 1} @current;
			$pare{$2} = [@current];
		}
	}
	
	#Distance matrix
	my %dis2;
	foreach my $i (0..$#order) {
		foreach my $j ($i..$#order) {
			if ($i == $j) {
				$dis2{$order[$i]}{$order[$j]} = 0;
			}
			else {
				my $tem = $dis{$order[$i]} + $dis{$order[$j]};
				my $tem2 = -1;
				foreach my $k (0..$#{$pare{$order[$i]}}) {
					last if ($k > $#{$pare{$order[$j]}});
					if ($pare{$order[$i]}[$k] eq $pare{$order[$j]}[$k]) {
						$tem2 = $k;
					}
				}
				
				if ($#{$pare{$order[$i]}} != -1) {
					map {$tem += $dis{'node_'.$_}} map {$pare{$order[$i]}[$_]} ($tem2+1)..$#{$pare{$order[$i]}};
				}
				if ($#{$pare{$order[$j]}} != -1) {
					map {$tem += $dis{'node_'.$_}} map {$pare{$order[$j]}[$_]} ($tem2+1)..$#{$pare{$order[$j]}};
				}
				$dis2{$order[$i]}{$order[$j]} = $dis2{$order[$j]}{$order[$i]} = $tem;
			}
		}
	}
	
	##output
	open (OUT, ">$out_file") or die "cannot open $out_file.\n";
	print OUT join("\t",'',@order),"\n";
	foreach my $i (@order) {
		print OUT join("\t",$i,map {$dis2{$i}{$_}} @order),"\n";
	}
	close(OUT);
}

#Given a FASTA file, generates a hash with all the sequences
sub fasta2hash {
	if (@_ != 1) {
		die "USAGE: module fasta2hash requires one arguments: <\$fasta_file>.\n";
	}
	
	my ($fasta_file) = @_;
	
	my %hash;
	open (FASTA, "<$fasta_file") or die "cannot open $fasta_file.\n";
	while (my $line = <FASTA>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			$header =~ s/^>//;
			
			my $seq = <FASTA>;
			$seq =~ s/\r//sig;
			$seq =~ s/\n//sig;
			
			$hash{$header} = $seq;
		}
		else {
			next;
		}
	}
	close (FASTA);
	
	my $hash = \%hash;
	return ($hash);
}

#Reads a list of headers and generates a FASTA file according to the database file.
sub fasta_maker {
	if (@_ != 3) {
		die "USAGE: module fasta_maker requires three arguments: <\@headers> <\%database> <\$out_file>.\n";
	}
	
	my ($headers, $database, $out_file) = @_;
	
	my @headers = @{ $headers };
	my %database = %{ $database };
		
	open (OUT_FILE, ">$out_file") or die "cannot open $out_file.\n";
	foreach my $header (sort @headers) {
		my $new_header = $header;
		$new_header =~ s/^>//;
		
		if (!exists $database{$new_header}) {
			die "ERROR: There is a header in the query that does not exist in the database.\n";
		}
		
		print OUT_FILE ">$new_header\n";
		print OUT_FILE "$database{$new_header}\n";
	}
	close (OUT_FILE);
}

#Reads a single-lined FASTA file and generates a FASTA file in a given directory for each sequence and returns the list of file names.
sub fasta_splitter {
	if (@_ != 2) {
		die "USAGE: module fasta_splitter requires two arguments: <\$in_file> <\$out_folder>.\n";
	}
	
	my ($in_file, $out_folder) = @_;
	my @out_files;
	
	open (IN, "<$in_file") or die "cannot open $in_file.\n";
	while (my $line = <IN>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			$header =~ s/^>//;
			my $seq = <IN>;
			$seq =~ s/\r//sig;
			$seq =~ s/\n//sig;
			
			#file name
			my $fname = $out_folder."/".$header.".fas";
			push (@out_files, $fname);
			#generating a file
			open (OUT, ">$fname") or die "cannot open $fname.\n";
			print OUT ">$header\n";
			print OUT "$seq\n";
			close (OUT);
		}
		else {
			next;
		}
	}
	close (IN);
	my $out_files = \@out_files;
	return ($out_files);
}

#given a tree, gene, and a threshold, cuts one branch that is longer than the threshold and gives back the subtree with the gene in it.
sub long_branch_cutter {
	if (@_ != 3) {
		die "USAGE: module long_branch_cutter requires three arguments: <\$test_tree>, <\$gene>, <\$dist_threshold>.\n";
	}
	
	my ($test_tree, $gene, $dist_threshold) = @_;
	
	#setting file names
	my $temp_tree_01 = $gene."_lbc01_temp.tre";
	my $temp_tree_02 = $gene."_lbc02_temp.tre";
	
	#make a file with the test_tree
	open (TEST_TREE, ">$temp_tree_01") or die "cannot open $temp_tree_01.\n";
	print TEST_TREE "$test_tree\;\n";
	close (TEST_TREE);
	
	#make a hash with all clades
	my %clades;
	my @list;
	push (@list, $test_tree);
	foreach my $element (@list) {
		my $mini_tree = $element;
		my @splitter = split (/\:/, $mini_tree);
		my $branch = $splitter[-1];
		if ($branch == 0) {
			$branch = 0.000000001;
		}
		pop (@splitter);
		$mini_tree = join ("\:", @splitter);
		$mini_tree =~ s/^\(//;
		$mini_tree =~ s/\)$//;
		$clades{$mini_tree} = $branch;
		
		if ($mini_tree =~ m/\,/) {
			my @mini_tree = split (//, $mini_tree);
			my $bracscore = 0;
			foreach my $word (@mini_tree) {
				if ($word =~ m/\(/) {
					$bracscore++;
					next;
				}
				elsif ($word =~ m/\)/) {
					--$bracscore;
					next;
				}
				elsif ($word =~ m/\,/ && $bracscore == 0) {
					$word = " ";
					next;
				}
				else {
					next;
				}
			}
			$mini_tree = join ("", @mini_tree);
			
			my @chunks = split (/ /, $mini_tree);
			foreach my $chunk (@chunks) {
				push (@list, $chunk);
			}
		}
		else {
			next;
		}
	}
	
	my $cut_test = 0;
	foreach my $clade (keys %clades) {
		if ($clades{$clade} > $dist_threshold) {
			#update test_tree
			if ($clade !~ m/$gene/) {
				open (TEMP, ">$temp_tree_02") or die "cannot open $temp_tree_02.\n";
				print TEMP "($clade)\:0\;\n";
				close (TEMP);
				my $removes = `nw_labels $temp_tree_02`;
				unlink ($temp_tree_02);
				$removes =~ s/\n$//;
				$removes =~ s/\r$//;
				my @removes = split (/\n/, $removes);
				$removes = join (" ", @removes);
				
				$test_tree = `nw_prune $temp_tree_01 $removes`;
				$test_tree =~ s/\r//sig;
				$test_tree =~ s/\n//sig;
				$test_tree =~ s/\;$//;
			}
			else {
				$test_tree = "(".$clade."):0";
			}
			$cut_test = 1;
			last;
		}
		else {
			next;
		}
	}
	
	unlink ($temp_tree_01);
	return ($test_tree, $cut_test);
}

#Reads a single-lined FASTA file and runs MAFFT to generate single-lined aligned file.
sub maffter {
	if (@_ != 5) {
		die "USAGE: module maffter requires five arguments: <\$sequence_file> <\$mafft_file> <\$sequence_type> <\$threads> <\$id>.\n";
	}
	
	#
	#Checking if MAFFT is running.
	#
	my $mafft_check = `which mafft`;
	if (!$mafft_check) {
		die "USAGE: module maffter requires MAFFT to be running.\n";
	}
	
	my ($sequence_file, $mafft_file, $sequence_type, $threads, $id) = @_;
	
	#setting parameters for mafft
	my $seqtype;
	if ($sequence_type =~ m/DNA/) {
		$seqtype = "--nuc";
	}
	elsif ($sequence_type =~ m/AA/) {
		$seqtype = "--amino";
	}
	else {
		die "ERROR: Sequence type needs to be DNA or AA.\n";
	}
	
	#testing for alignment methods
	my $seq_count = `cat $sequence_file | wc -l`;
	$seq_count =~ s/\r//sig;
	$seq_count =~ s/\n//sig;
	my @seq_count = split (" ", $seq_count);
	$seq_count = $seq_count[-1] / 2;
	
	#mafft conditions
	my $max_iteration = 1000;
	my $retree = 2;
	my $aln_method;
	if ($seq_count < 200) {
		$aln_method = "--localpair";
	}
	else {
		$aln_method = "--6merpair";
	}
	
	#temp files
	my $temp1 = "MAFFTER_01_".$id.".fas";
	my $temp2 = "MAFFTER_02_".$id.".fas";
	
	#running MAFFT
	system "mafft $seqtype $aln_method --retree $retree --maxiterate $max_iteration --anysymbol --quiet --thread $threads $sequence_file > $temp1";
	
	&OrthoReD_library::singler ($temp1, $temp2);
	
	open (PRE, "<$temp2") or die "cannot open $temp2.\n";
	open (POST, ">$mafft_file") or die "cannot open $mafft_file.\n";
	while (my $line = <PRE>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		if ($line =~ m/^>/) {
			print POST "$line\n";
			next;
		}
		else {
			$line =~ s/U/X/g;
			print POST "$line\n";
		}
	}
	close (PRE);
	close (POST);
	
	unlink ($temp1);
	unlink ($temp2);
	
	my $temp_log;
	$temp_log = "Stop codons\(\*\) are not included in the alignement.\n";
	$temp_log = $temp_log."Selenocysteine(U) is replaced by X.\n";
	
	return ($temp_log);
}

#Given a BLAST hit table, generates a file set up for running MCL
sub mcl_prepper {
	if (@_ != 3) {
		die "USAGE: module mcl_prepper requires three arguments: <\$in_file> <\$out_file> <\$id>.\n";
	}
	
	my ($in_file, $out_file, $id) = @_;
	
	#prep_file
	my $temp_01 = "MCL_PREPPER_01_".$id.".txt";
	
	system "cut -f 1-3 $in_file > $temp_01";
	
	#generating a list of unique members of the BLAST table
	my $members = `cut -f 1 $temp_01`.`cut -f 2 $temp_01`;
	$members =~ s/\n$//;
	$members =~ s/\r$//;
	my @members = split (/\n/, $members);
	
	my %seen;
	my @unique;
	foreach my $member (@members) {
		if (! $seen{$member}) {
			push (@unique, $member);
			$seen{$member} = 1;
		}
	}
	@members = sort @unique;
	
	#generating a two-dimensional array of members with default values
	my @matrix;
	my %members;
	my $order = 0;
	foreach my $member (@members) {
		$members{$member} = $order;
		$order++;
	}
	for (my $i=0; $i < $order; $i++) {
		for (my $j=0; $j < $order; $j++) {
			$matrix[$i][$j] = "9.0e-01";
		}
	}
	
	#adding actual values into the matrix
	open (PREP, "<$temp_01") or die "cannot open $temp_01.\n";
	while (my $line = <PREP>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		my @line = split (/\t/, $line);
		my @sq = ($line[0], $line[1]);
		@sq = sort @sq;
		
		if ($line[2] < $matrix[$members{$sq[0]}][$members{$sq[1]}]) {
			$matrix[$members{$sq[0]}][$members{$sq[1]}] = $line[2];
		}
		else {
			next;
		}
	}
	close (PREP);
	
	#printing out a file prepared for MCL
	open (OUT, ">$out_file") or die "cannot open $out_file.\n";
	for (my $i=0; $i < $order; $i++) {
		for (my $j=0; $j < $order; $j++) {
			my @unsort = ($members[$i], $members[$j]);
			my $unsort = join ("", @unsort);
			my @sort = sort @unsort;
			my $sort = join ("", @sort);
			
			if ($sort !~ m/$unsort/) {
				next;
			}
			else {
				print OUT "$members[$i]\t$members[$j]\t$matrix[$i][$j]\n";
			}
		}
	}
	close (OUT);
	
	unlink ($temp_01);
}

#Given a distance matrix and parameters, runs mcl and returns a list of all members of the cluster with the query.
sub mcl_runner {
	if (@_ != 5) {
		die "USAGE: module mcl_runner requires four arguments: <\$matrix_file> <\$query> <\$param> <\$id> <\$threads>.\n";
	}
	
	my ($matrix_file, $query, $param, $id, $threads) = @_;
	
	#out_file
	my $temp_01 = "MCL_RUNNER_01_".$id.".txt";
	my $temp_02 = "MCL_RUNNER_02_".$id.".mci";
	my $temp_03 = "MCL_RUNNER_03_".$id.".tab";
	

	
	my $mcl_temp = "temp_mcl.txt";
	system "mcxload -abc $matrix_file --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $temp_02 -write-tab $temp_03 2>/dev/null";
	system "mcl $temp_02 -I $param -use-tab $temp_03 -o $temp_01 -q x -V all -te $threads 2>/dev/null";
	
	my $cluster = `grep $query $temp_01`;
	$cluster =~ s/\n/\t/g;
	$cluster =~ s/\r/\t/g;
	$cluster =~ s/\t$//g;
	my @cluster = split(/\t/, $cluster);
	$cluster = \@cluster;
	
	unlink ($temp_01);
	unlink ($temp_02);
	unlink ($temp_03);
	
	return ($cluster);

}

#Given a input tree file and an output tree file, roots the input tree at midpoint and generates the output tree.
sub midpointer {
	if (@_ != 2) {
		die "USAGE: module midpointer requires three arguments: <\$in_tree>, <\$out_tree>.\n";
	}
	
	my ($in_tree, $out_tree) = @_;
	my $dist_file = $in_tree;
	$dist_file =~ s/.tre$/_MIDPOINT_DIST.txt/;
	
	#making a distance matrix of the in_tree
	&distance_finder ($in_tree, $dist_file);
	
	#reading in the tip names in the matrix
	my $tips = `head -n 1 $dist_file`;
	$tips =~ s/^\t//;
	$tips =~ s/\r//sig;
	$tips =~ s/\n//sig;
	my @tips = split (/\t/, $tips);
	my %tips;
	my $order = 1;
	for my $element (@tips) {
		$tips{$order} = $element;
		$order++;
	}
	
	#identifying the most distant tips and their distance
	my @dist_pair;
	my $max_dist = -1;
	my $threshold_dist;
	
	my $line_count = 0;
	open (DIST, "<$dist_file") or die "cannot open $dist_file.\n";
	while (my $line = <DIST>) {
		$line_count++;
		if ($line_count > 1) {
			my $dist = $line;
			$dist =~ s/\r//sig;
			$dist =~ s/\n//sig;
			my @dist = split (/\t/, $dist);
			shift @dist;
			
			my $dist_count = 0;
			foreach my $element (@dist) {
				$dist_count++;
				if ($element > $max_dist) {
					$max_dist = $element;
					
					@dist_pair = ();
					push (@dist_pair, $tips{$line_count - 1});
					push (@dist_pair, $tips{$dist_count});
				}
				else {
					next;
				}
			}
		}
		else {
			next;
		}
	}
	close (DIST);
	$threshold_dist = $max_dist/2;
	unlink ($dist_file);
	
	#in case the brach lengths are all 0.
	if ($threshold_dist == 0) {
		system "cp $in_tree $out_tree";
	}
	else {
		#root the tree temporarily based on a tip that is not the most distant pair
		my $temproot1 = `nw_labels $in_tree | head -n 3`;
		$temproot1 =~ s/\r$//;
		$temproot1 =~ s/\n$//;
		
		my @temproot1 = split (/\n/, $temproot1);
		if ($temproot1[0] !~ m/$dist_pair[0]/ && $temproot1[0] !~ m/$dist_pair[1]/) {
			$temproot1 = $temproot1[0];
		}
		elsif ($temproot1[1] !~ m/$dist_pair[0]/ && $temproot1[1] !~ m/$dist_pair[1]/) {
			$temproot1 = $temproot1[1];
		}
		elsif ($temproot1[2] !~ m/$dist_pair[0]/ && $temproot1[2] !~ m/$dist_pair[1]/) {
			$temproot1 = $temproot1[2];
		}
		else {
			die "ERROR: temporary root is not properly selected.\n";
		}
		
		my $temproot_tree1 = `nw_reroot $in_tree $temproot1`;
		$temproot_tree1 =~ s/\r//sig;
		$temproot_tree1 =~ s/\n//sig;
		
		#root it according to one of the most distant pair with the shorter branch length
		my $pair1_bl = `echo '$temproot_tree1' | nw_distance -m p - $dist_pair[0] | cut -f 2`;
		$pair1_bl =~ s/\r//sig;
		$pair1_bl =~ s/\n//sig;
		
		my $pair2_bl = `echo '$temproot_tree1' | nw_distance -m p - $dist_pair[1] | cut -f 2`;
		$pair2_bl =~ s/\r//sig;
		$pair2_bl =~ s/\n//sig;
		
		if ($pair2_bl < $pair1_bl) {
			@dist_pair = reverse @dist_pair;
		}
		
		my $temproot_tree2 = `echo '$temproot_tree1' | nw_reroot - $dist_pair[0]`;
		$temproot_tree2 =~ s/\r//sig;
		$temproot_tree2 =~ s/\n//sig;
		
		#find the clade that would be used for midpoint rooting
		my $test_tre = $temproot_tree2;
		my $cycle = 0;
		my $branch_total = 0;
		my @keep_tips = @tips;
		while ($branch_total <= $threshold_dist) {
			$cycle++;
			
			$test_tre =~ s/\;$//;
			$test_tre =~ s/^\(//;
			$test_tre =~ s/\)$//;
			
			my @test_tre;
			if ($test_tre =~ m/\,/) {
				@test_tre = split (//, $test_tre);
				my $bracscore = 0;
				foreach my $word (@test_tre) {
					if ($word =~ m/\(/) {
						$bracscore++;
						next;
					}
					elsif ($word =~ m/\)/) {
						--$bracscore;
						next;
					}
					elsif ($word =~ m/\,/ && $bracscore == 0) {
						$word = " ";
						next;
					}
					else {
						next;
					}
				}
			}
			else {
				push (@test_tre, $test_tre);
			}
			$test_tre = join ("", @test_tre);
			
			my @chunks = split (/ /, $test_tre);
			
			foreach my $chunk (@chunks) {
				if ($chunk =~ m/$dist_pair[1]/) {
					$test_tre = $chunk;
					last;
				}
				else {
					next;
				}
			}
			
			my @splitter = split (/\:/, $test_tre);
			my $branch = $splitter[-1];
			pop (@splitter);
			if ($cycle == 1) {
				$branch = $branch * 2;
			}
			$branch_total = $branch_total + $branch;
			
			$test_tre = join ("\:", @splitter);
			$test_tre = $test_tre."\;";
			
			my $keep_tips = `echo '$test_tre' | nw_labels -It -`;
			$keep_tips =~ s/\r//sig;
			$keep_tips =~ s/\n//sig;
			@keep_tips = split (/\t/, $keep_tips);
		}
		
		#root the tree according to the determined tips
		my $midpoint_tips = join (" ", @keep_tips);
		my $midpoint_tre = `nw_reroot -l $in_tree $midpoint_tips`;
		if (!$midpoint_tre) {
			$midpoint_tre = $in_tree;
		}
		$midpoint_tre =~ s/\r//sig;
		$midpoint_tre =~ s/\n//sig;
		
		open (OUT, ">$out_tree") or die "cannot open $out_tree.\n";
		print OUT "$midpoint_tre\n";
		close (OUT);
	}
}

#Reads a single-lined FASTA file with three of less sequences and generates a parenthetical tree in NEXUS format.
sub mini_tree_maker {
	if (@_ != 2) {
		die "USAGE: module mini_tree_maker requires two arguments: <\$sequence_file> <\$tree_file>.\n";
	}
	
	my ($sequence_file, $tree_file) = @_;
	
	my @tips;
	open (FASTA, "<$sequence_file") or die "cannot open $sequence_file.\n";
	while (my $line = <FASTA>) {
		if ($line =~ m/^>/) {
			$line =~ s/\r//sig;
			$line =~ s/\n//sig;
			my $header = $line;
			$header =~ s/^>//;
			push(@tips, $header);
		}
	}
	close (FASTA);
	my $num = @tips;
	
	open (TREE, ">$tree_file") or die "cannot open $tree_file.\n";
	if ($num == 1) {
		print TREE "$tips[0]\;\n";
	}
	elsif ($num == 2) {
		print TREE "($tips[0]:1,$tips[1]:1)\;\n";
	}
	elsif ($num == 3) {
		print TREE "($tips[0]:1,$tips[1]:1,$tips[2]:1)\;\n";
	}
	else {
		die "ERROR: mini_tree_maker cannot handle tree with more than three tips.\n";
	}
	close (TREE);
}

#Given two arrays, generate tree sets of arrays; arrays with unique elements in each array, and an array with all the common elements.
sub mismatch_finder {
	if (@_ != 2) {
		die "USAGE: module mismatch_finder requires two arguments: <\$array1> <\$array2>.\n";
	}
	
	my @array1 = @{ $_[0] };
	my @array2 = @{ $_[1] };
	
	#Making an array of "all unique elements", "common", and "difference".
	my @all;
	my @common;
	my @difference;
	my @only1;
	my @only2;
	
	my %count1;
	foreach my $element (@array1, @array2) {
		push (@all, $element);
		$count1{$element}++;
	}
	foreach my $element (keys %count1) {
		if ($count1{$element} > 1) {
			push (@common, $element);
			next;
		}
		push (@difference, $element);
	}
	
	my %count2;
	foreach my $element (@array1, @difference) {
		$count2{$element}++;
	}
	foreach my $element (keys %count2) {
		if ($count2{$element} > 1) {
			push (@only1, $element);
		}
	}
	
	my %count3;
	foreach my $element (@array2, @difference) {
		$count3{$element}++;
	}
	foreach my $element (keys %count3) {
		if ($count3{$element} > 1) {
			push (@only2, $element);
		}
	}
	
	my $common = \@common;
	my $only1 = \@only1;
	my $only2 = \@only2;
	
	return ($common, $only1, $only2);
}

#Reads a FASTA file with DNA sequences to check if it is in single-lined FASTA format with full length DCS.
sub qcer_file {
	if (@_ != 3) {
		die "USAGE: module quer_file requires three argument: <\$in_file> <\$sequence_type> <\$full_length>.\n";
	}
	
	my ($in_file, $sequence_type, $full_length) = @_;
	my $value = "TRUE";
	
	open (QC, "<$in_file") or die "cannot open $in_file.\n";
	my $line_count = 0;
	my $header;
	my $sequence;
	while (my $line = <QC>) {
		$line_count++;
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		
		if ($line =~ m/^$/) {
				$value = "ERROR: $in_file has a blank line. Check line $line_count.";
				last;
		}
		
		if ($line_count % 2 == 1) {
			if ($line !~ m/^>/) {
				$value = "ERROR: $in_file is not in single-line FASTA format. Check line $line_count.";
				last;
			}
			$header = $line;
			($value) = &OrthoReD_library::qcer_line ("header", $header, $sequence_type, $full_length);
			if ($value !~ m/^TRUE$/) {
				$value = $value." Check line $line_count.";
				last;
			}
		}
		
		if ($line_count % 2 == 0) {
			$sequence = $line;
			if ($line =~ m/^>/) {
				$value = "ERROR: $in_file is not in single-line FASTA format. Check line $line_count.";
				last;
			}
			($value) = &OrthoReD_library::qcer_line ("sequence", $sequence, $sequence_type, $full_length);
			if ($value !~ m/^TRUE$/) {
				$value = $value." Check line $line_count.";
				last;
			}
		}
	}
	close (QC);
	return ($value);
}

#Reads the header or the sequence of a DNA FASTA file and checks ifs quality.
	#Any header that does not begin with ">" or includes "-", ":", or "|" will be reported.
	#Any sequences that are not full-length coding sequence will be reported.
		#Sequences need to begin with ATG.
		#Sequences need to be 3n long.
		#Sequences cannot include gaps (-).
		#Sequences cannot include stop codons other than at the end.
		#Sequence needs to end with a stop codon.
sub qcer_line {
	if (@_ != 4) {
		die "USAGE: module quality_checker requires four arguments: <header/sequence> <\$pre_qc> <\$sequence_type> <\$full_length>.\n";
	}
	
	my ($line_type) = $_[0];
	my ($header) = $_[1];
	my ($sequence) = $_[1];
	my ($sequence_type) = $_[2];
	my ($full_length) = $_[3];
	
	my %translator1 = &OrthoReD_library::translator1 ();
	my $value = "TRUE";
	
	#
	#Handling headers
	#
	if ($line_type =~ m/^header$/) {
		$header =~ s/\r//sig;
		$header =~ s/\n//sig;
		
		if ($header !~ m/^>/) {
			$value = "ERROR: Header not in single-lined FASTA format.";
		}
		elsif ($header =~ m/-/) {
			$value = "ERROR: Headers cannot include the symble (-).";
		}
		elsif ($header =~ m/\:/) {
			$value = "ERROR: Headers cannot include the symble (:).";
		}
		elsif ($header =~ m/\|/) {
			$value = "ERROR: Headers cannot include the symble (|).";
		}
		elsif ($header =~ m/ /) {
			$value = "ERROR: Headers cannot include the symble ( ).";
		}
	}
	
	#
	#Handling sequence
	#
	elsif ($line_type =~ m/^sequence$/ and $sequence_type =~ m/DNA/) {
		$sequence =~ s/\r//sig;
		$sequence =~ s/\n//sig;
		
		if ($sequence =~ m/>/) {
			$value = "ERROR: Sequence not in single-lined FASTA format.";
		}
		my $sequence_length = length ($sequence);
		if ($sequence_length % 3 != 0) {
			$value = "ERROR: Sequence needs to be 3n long.";
		}
		if ($sequence =~ m/-/) {
			$value = "ERROR: Sequence cannot include gaps (-).";
		}
		if ($sequence !~ m/^ATG/ and $full_length =~ m/YES/) {
			$value = "ERROR: Sequence must start with ATG.";
		}
		for (my $c = 1; $c <= ($sequence_length - 3) / 3; $c++) {
			my $codon = substr ($sequence, $c * 3 - 3, 3);
			if (exists $translator1{$codon} and $translator1{$codon} =~ m/\*/) {
				$value = "ERROR: Sequence cannot include stop codons other than at the end.";
			}
		}
		if ($sequence !~ m/TAA$/ and $sequence !~ m/TAG$/ and $sequence !~ m/TAR$/ and $sequence !~ m/TGA$/ and $sequence !~ m/TRA$/ and $full_length =~ m/YES/) {
			$value = "ERROR: Sequence must end with a stop codon.\n";
		}
	}
	elsif ($line_type =~ m/^sequence$/ and $sequence_type =~ m/AA/) {
		$sequence =~ s/\r//sig;
		$sequence =~ s/\n//sig;
		
		if ($sequence =~ m/^>/) {
			$value = "ERROR: Sequence looks like a header starting with a (\>).";
		}
		
		my $tempseq = $sequence;
		$tempseq =~ s/\*$//g;
		
		if ($tempseq =~ m/\*/) {
			$value = "ERROR: Sequence cannot include character (\*) other than at the end of a sequence.";
		}
		elsif ($tempseq =~ m/[^a-zA-Z]/){
			$value = "ERROR: Sequence cannot include non-alphabet characters other than (\*).";
		}
		elsif ($tempseq =~ m/[BJOUZ]/i) {
			$value = "ERROR: Sequence contains characters that are not recognized as regular amino acids (B, J, O, U, or Z).";
		}
		elsif ($tempseq !~ m/^M/ and $full_length =~ m/YES/) {
			$value = "ERROR: Sequence must start with M.";
		}
	}
	
	else {
		$value = "ERROR: Invalid options were given to module qcer_line.";
	}
	
	return ($value);
}

#Given a BLAST output file and a directory, splits the BLAST output into individual queries
sub query_separator {
	if (@_ != 2) {
		die "USAGE: module query_separator requires two arguments: <\$blastfile>, <\$folder>.\n";
	}
	
	my ($blastfile, $folder) = @_;
	
	my $inquerys = `cut -f1 $blastfile | sort | uniq`;
	$inquerys =~ s/\r$//;
	$inquerys =~ s/\n$//;
	my @inquerys = split (/\n/, $inquerys);
	@inquerys = sort @inquerys;
	
	foreach my $inquery (@inquerys) {
		my $outfile = $folder."/".$inquery;
		system "grep '^$inquery' $blastfile >$outfile";
	}
	
	$inquerys = "";
	$inquerys = \@inquerys;
	
	return ($inquerys);
}

#given two directories, one with all orthologs and another with BLAST hits, ranks each orthologs and summarises them.
sub ranker {
	if (@_ != 3) {
		die "USAGE: module ranker requires three arguments: <\$ortho_folder> <\$hit_folder> <\$out_file>.\n";
	}
	
	my ($ortho_folder, $hit_folder, $out_file) = @_;
	
	#making a hash for ranking
	my $total = 0;
	my %rankcount;
	
	#make a list of files to look at
	my $keys = `ls $ortho_folder`;
	$keys =~ s/\r$//;
	$keys =~ s/\n$//;
	my @keys = split (/\n/, $keys);
	@keys = sort @keys;
	
	#check the ranking of each key
	foreach my $key (@keys) {
		my @elements = split (/_/, $key);
		$key = $elements[0]."_".$elements[1]."_".$elements[2]."_".$elements[3];
		
		my $ortho_file = `ls $ortho_folder | grep '$key'`;
		$ortho_file = $ortho_folder."/".$ortho_file;
		$ortho_file =~ s/\r//sig;
		$ortho_file =~ s/\n//sig;
		my $ortholist = `grep '>' $ortho_file | sed 's/>//'`;
		$ortholist =~ s/\r$//;
		$ortholist =~ s/\n$//;
		my @ortholist = split (/\n/, $ortholist);
		
		my $hit_file = `ls $hit_folder | grep '$key'`;
		$hit_file = $hit_folder."/".$hit_file;
		$hit_file =~ s/\r//sig;
		$hit_file =~ s/\n//sig;
		
		my $uniq_hit_file1 = "RANKER_".$key."_01\.txt";
		
		&OrthoReD_library::uniq_finder($hit_file, $uniq_hit_file1);
		
		foreach my $ortholog (@ortholist) {
			$total++;
			my $sp = $ortholog;
			my @sp = split (/_/, $sp);
			$sp = $sp[0];
			my $locus = $sp[1];
			
			my $rank;
			my $eval = `cut -f2-3 $uniq_hit_file1 | grep $ortholog | cut -f2`;
			if (!$eval) {
				$rank = 1;
			}
			elsif ($eval =~ m/^0$/ or $eval =~ m/^0\.$/ or $eval =~ m/^0\.0$/) {
				$rank = 1;
			}
			else {
				my $cutoff = `cut -f2-3 $uniq_hit_file1 | grep ^$sp | sort -k2,2g | grep -n $locus | head -n 1 | cut -d '\:' -f1`;
				chomp $cutoff;
				my $value = `grep $locus $uniq_hit_file1 | head -n 1 | cut -f3`;
				chomp $value;
				my $tie = `cut -f2-3 $uniq_hit_file1 | grep ^$sp | head -n $cutoff | grep $value | wc -l`;
				$tie =~ s/ //g;
				chomp $tie;
				$rank = $cutoff - $tie + 1;
			}
			
			if (!exists $rankcount{$rank}) {
				$rankcount{$rank} = 1;
			}
			else {
				$rankcount{$rank}++;
			}
		}
		unlink ($uniq_hit_file1);
	}
	
	#print out the rankings
	open (OUT, ">$out_file") or die "cannot open $out_file.\n";
	print OUT "Total pridicted orthologs:\t$total\n";
	foreach my $key (sort {$a <=> $b} keys %rankcount) {
		my $percentage = sprintf ("%.2f", ($rankcount{$key} / $total) * 100);
		print OUT "No.$key ranked:\t$rankcount{$key}($percentage\%)\n";
	}
	close (OUT);
}

#Given a single-lined FASTA file of AA, runs RAXML to generate a tree with branchlength as quickly as possible
sub raxmler_speed {
	if (@_ != 5) {
		die "USAGE: module raxmler requires five arguments: <\$in_file> <\$out_file> <\$threads> <\$id> <\$vraxml>.\n";
	}
	
	my ($in_file, $out_file, $threads, $id, $vraxml) = @_;
	my $return = "GOOD";
	
	#determining the correct version of raxml.
	my $command;
	if ($vraxml =~ m/^AVX$/) {
		$command = "raxmlHPC-PTHREADS-AVX";
	}
	elsif ($vraxml =~ m/^SSE3$/) {
		$command = "raxmlHPC-PTHREADS-SSE3";
	}
	else {
		die "ERROR: module raxmler requires \$vraxml to be AVX or SSE3.\n";
	}
	
	#Specifying model for RAXML run
	my $m = "PROTCATIAUTO";
	my $aminotest = &OrthoReD_library::aminotester ($in_file);
	if ($aminotest =~ m/FAIL/) {
		$m = "PROTCATIJTTX";
	}
	
	#out_file
	my $out = "RAXMLER_SPEED_".$id;
	mkdir ($out);
	
	#make tree quickly
	system "$command -f E -F -s $in_file -m $m -p 12345 -n QUICK_$id -T $threads";
	system "mv *QUICK_$id* $out";
	
	#checking if backup tree making is needed.
	my $test = `[ -f $out/RAxML_fastTree.QUICK_$id ] && echo "Found" || echo "Not found"`;
	if ($test =~ m/Not found/) {
		$m = "PROTCATIJTTX";
		system "rm $out/*";
		system "$command -f E -F -s $in_file -m $m -p 12345 -n QUICK_$id -T $threads";
		system "mv *QUICK_$id* $out";
	}
	
	#checking if the backup worked
	$test = `[ -f $out/RAxML_fastTree.QUICK_$id ] && echo "Found" || echo "Not found"`;
	if ($test =~ m/Not found/) {
		my $tips = `grep '>' $in_file`;
		$tips =~ s/\r/\t/g;
		$tips =~ s/\n/\t/g;
		$tips =~ s/\t\t/\t/g;
		$tips =~ s/\>//g;
		$tips =~ s/\t/\:1\,/g;
		$tips =~ s/\,$//;
		$tips = "((".$tips.")\:0)\;";
		
		open(OUT, ">$out_file") or die "cannot open $out_file.\n";
		print OUT "$tips\n";
		close(OUT);
		$return = "BACKUP";
	}
	else {
		#adding branchlength to the tree
		system "$command -f e -t $out/RAxML_fastTree.QUICK_$id -m $m -s $in_file -n QUICK_BL_$id -T $threads";
		system "mv *QUICK_BL_$id* $out";
		system "cp $out/RAxML_result.QUICK_BL_$id $out_file";
	}
	system "rm -rf $out";
	return ($return);
}

#Given a tree and a root, roots the tree at the root (or at midpoint)
sub rooter {
	if (@_ != 3) {
		die "USAGE: module rooter requires three arguments: <\$tree_u> <\$root> <\$tree_r>.\n";
	}
	
	my ($tree_u, $root, $tree_r) = @_;
	
	#testing if root exists in the tree
	my $tre_line = `head -n 1 $tree_u`;
	$tre_line =~ s/\r//sig;
	$tre_line =~ s/\n//sig;
	
	if ($tre_line !~ m/$root/ and $root !~ m/^midpoint$/) {
		die "ERROR: $root is not a tip present in $tree_u.\n";
	}
	
	if ($root =~ m/^midpoint$/) {
		&midpointer ($tree_u, $tree_r);
	}
	else {
		system "nw_reroot $tree_u $root > $tree_r";
	}
}

#Given a gene of interest, distance matrix of the tree and a list of outgroups, picks on tip as its root. (midpoint when outgroup is absent)
sub root_finder {
	if (@_ != 6) {
		die "USAGE: module root_finder requires six arguments: <\$goi> <\$tree> <\$dist> <\$og> <\$spp_list> <\$rooting>.\n";
	}
	
	my ($goi, $tree, $dist, $ogs, $spp_list, $rooting) = @_;
	
	my $full_name = 0;
	my $short_name = 1;
	my $acronym_name = 2;
	my %spp = %{ &OrthoReD_library::species_lister ($spp_list, $full_name, $acronym_name) };
	
	my @ogs = @{ $ogs };
	foreach my $og (@ogs) {
		$og = $spp{$og};
	}
	my %ogs = map { $_ => 1 } @ogs;
	
	#selecting the outgroup tips
	my @ogtips;
	foreach my $og (@ogs) {
		my $tips = `nw_labels $tree | grep $og`;
		$tips =~ s/\r$//;
		$tips =~ s/\n$//;
		my @tips = split (/\n/, $tips);
		push (@ogtips, @tips);
	}
	my $ogtips = scalar (@ogtips);
	
	#rooting
	my $root = "midpoint";
	#if no outgroup is available
	if ($ogtips < 1) {
		return ($root);
	}
	else {
		#making a distance matrix of the trees
		&OrthoReD_library::distance_finder ($tree, $dist);
		
		#making a list of tips with their IDs
		my %tip_id;
		my $tips = `head -n 1 $dist`;
		$tips =~ s/\r//sig;
		$tips =~ s/\n//sig;
		$tips =~ s/^\t//;
		my @tips = split (/\t/, $tips);
		
		#selecting the line with the gene of interest;
		my $goi_line = `grep ^$goi $dist | tail -1`;
		$goi_line =~ s/\r//sig;
		$goi_line =~ s/\n//sig;
		my @goi_dist = split (/\t/, $goi_line);
		
		#linking the tips and distances
		my %og_dist;
		foreach my $og (@ogtips) {
			$og_dist{$og} = 0;
		}
		
		my $element_count = 1;
		foreach my $tip (@tips) {
			if (exists $og_dist{$tip}) {
				$og_dist{$tip} = $goi_dist[$element_count];
			}
			$element_count++;
		}
		
		#if rooting option is 'MD', most distant outgroup
		if ($rooting =~ m/MD/) {
			my $root_dist = 0;
			foreach my $tip (sort keys %og_dist) {
				if ($root_dist < $og_dist{$tip}) {
					$root = $tip;
					$root_dist = $og_dist{$tip};
				}
				else {
					next;
				}
			}
			return ($root);
		}
		#if rooting option is 'MI', most inclusve, most distant outgroup
		elsif ($rooting =~ m/MI/) {
			my $max_inspp = 0;
			my %inspp;
			foreach my $tip (sort keys %og_dist) {
				my $test = 0;
				my $upnode = 0;
				while ($test < 1) {
					my $spp = `nw_reroot $tree $tip | nw_clade -c $upnode - $goi | nw_labels - | cut -f1 -d'_' | sort | uniq`;
					$spp =~ s/\r$//;
					$spp =~ s/\n$//;
					my @spp = split (/\n/, $spp);
					$upnode++;
					
					foreach my $sp (@spp) {
						if (exists $ogs{$sp}) {
							$test++;
						}
					}
					
					if ($test < 1) {
						my @inspp = @spp;
						my $inspp = scalar (@inspp);
						$inspp{$tip} = $inspp;
					}
				}
			}
			
			foreach my $tip (sort keys %inspp) {
				if ($inspp{$tip} > $max_inspp) {
					$root = $tip;
					$max_inspp = $inspp{$tip};
				}
				elsif ($inspp{$tip} == $max_inspp && $og_dist{$tip} > $og_dist{$root}) {
					$root = $tip;
					$max_inspp = $inspp{$tip};
				}
				else {
					next;
				}
			}
			return ($root);
		}
	}
}

#Reads a single-lined FASTA file of DNA and generates two single-lined FASTA files: one with only full-length CDS, another with all removed sequences.
#Additionally it provides the % sequence removed.
#Any sequences that are not full-length coding sequence will be removed.
	#Sequences need to begin with ATG.
	#Sequences need to be 3n long.
	#Sequences cannot include gaps (-).
	#Sequences cannot include stop codons other than at the end.
	#Sequence needs to end with a stop codon.
sub separator {
	if (@_ != 6) {
		die "USAGE: module separater requires six arguments: <\$species_list> <\$in_file> <\$cds_file> <\$removed_file> <\$sequence_tpye> <\$full_length>.\n";
	}
	
	my ($species_list, $file_name, $cds_file, $removed_file, $sequence_type, $full_length) = @_;
	
	my $full_name = 0;
	my $short_name = 1;
	my $acronym_name = 2;
	my %spp = %{ &OrthoReD_library::species_lister ($species_list, $acronym_name, $full_name) };
	
	my %translator1 = &OrthoReD_library::translator1 ();
	
	my $percent_removed;
	my $total_line_count = 0;
	my $removed_line_count = 0;
	
	open (SINGLE, "<$file_name") or die "cannot open $file_name.\n";
	open (CDS, ">$cds_file") or die "cannot open $cds_file.\n";
	open (REMOVED, ">$removed_file") or die "cannot open $removed_file.\n";
	while ($_ = <SINGLE>) {
		if ($_ =~ m/^>/) {
			my $header = $_;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			
			my $cds = <SINGLE>;
			$cds =~ s/\r//sig;
			$cds =~ s/\n//sig;
			
			my ($header_value) = &OrthoReD_library::qcer_line ("header", $header, $sequence_type, $full_length);
			my ($seq_value) = &OrthoReD_library::qcer_line ("sequence", $cds, $sequence_type, $full_length);
			if ($header_value !~ m/^TRUE$/ or $seq_value !~ m/^TRUE$/) {
				print REMOVED "$header\n";
				print REMOVED "$cds\n";
				$total_line_count++;
				$removed_line_count++;
				next;
			}
			print CDS "$header\n";
			print CDS "$cds\n";
			$total_line_count++;
		}
		else {
			next;
		}
	}
	close (SINGLE);
	close (CDS);
	close (REMOVED);
	$percent_removed = sprintf ("%.2f", ($removed_line_count / $total_line_count) * 100);
	return ($percent_removed);
}

#Given an alignment file, removes identical sequences and generates a reduced file.
sub seq_reducer {
	if (@_ != 3) {
		die "USAGE: module seq_reducer requires three arguments: <\$full_file> <\$short_file> <\$conv_list>.\n";
	}
	
	my ($full_file, $short_file, $conv_list) = @_;
	
	#Reading the sequence file
	my %sequences;
	
	open (SEQS, "<$full_file") or die "cannot open $full_file.\n";
	while (my $line = <SEQS>){
		$line =~ s/\n//sig;
		$line =~ s/\r//sig;
		if ($line =~ m/^>/) {
			my $header = $line;
			$header =~ s/^>//;
			my $sequence = <SEQS>;
			$sequence =~ s/\n//sig;
			$sequence =~ s/\r//sig;
			
			my @headers = ();
			if (!exists $sequences{$sequence}) {
				$headers[0] = $header;
				$sequences{$sequence} = \@headers;
				next;
			}
			else {
				@headers = @{ $sequences{$sequence} };
				push (@headers, $header);
				$sequences{$sequence} = \@headers;
				next;
			}
		}
	}
	close (SEQS);
	
	#print out the conversion list
	open (LIST, ">$conv_list") or die "cannot open $conv_list.\n";
	foreach my $sequence (sort keys %sequences) {
		my @headers = @{ $sequences{$sequence} };
		my $headers = join ("\t", @headers);
		print LIST "$headers\n";
		
		$sequences{$sequence} = $headers[0];
	}
	close (LIST);
	
	#print out unique sequences.
	my %forprint;
	foreach my $sequence (sort keys %sequences) {
		$forprint{ $sequences{$sequence} } = $sequence;
	}
	
	open (UNQ, ">$short_file") or die "cannot open $short_file.\n";
	foreach my $header (sort keys %forprint) {
		print UNQ ">$header\n";
		print UNQ "$forprint{$header}\n";
	}
	close (UNQ);
}

#Reads a multi-lined FASTA file and generates a single-lined FASTA file.
sub singler {
	if (@_ != 2) {
		die "USAGE: module singler requires two arguments: <\$multi_file> <\$single_file>.\n";
	}
	
	my ($multi_file, $single_file) = @_;
	
	open (MULTI, "<$multi_file") or die "cannot open $multi_file.\n";
	open (SINGLE, ">$single_file") or die "cannot open $single_file.\n";
	my $line_count = 0;
	while (my $line = <MULTI>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		$line_count++;
		if ($line_count == 1 and $line =~ m/^>/) {
			print SINGLE "$line\n";
			next;
		}
		
		if ($line =~ m/^>/) {
			print SINGLE "\n$line\n";
			next;
		}
		else {
			$line = uc($line);
			print SINGLE "$line";
		}
	}
	print SINGLE "\n";
	close (MULTI);
	close (SINGLE);
}

#Reads a species list and generates a hash linking one column with another.
sub species_lister {
	if (@_ != 3) {
		die "USAGE: module species_lsiter requires three arguments:<\$species_list> <\$key> <\$value>\n";
	}
	
	my ($species_list, $key, $value) = @_;
	
	my %hash;
	open (SPP, "<$species_list") or die "cannot open $species_list.\n";
	while (my $sp = <SPP>) {
		$sp =~ s/\r//sig;
		$sp =~ s/\n//sig;
		my @line = split (/\t/, $sp);
		$hash{$line[$key]} = $line[$value];
	}
	close (SPP);
	return (\%hash);
}

#Reads set of FASTA files and generates a table summarizing them.
sub summarizer {
	if (@_ != 3 and @_ != 4) {
		die "USAGE: module summarizer requires three or four arguments: <\$species_list> <\$in_folder> <\$out_file> <\$conversion_list(optional)>.\n";
	}
	
	my ($species_list) = $_[0];
	my ($in_folder) = $_[1];
	my ($out_file) = $_[2];
	my $conversion_list = "NONE";
	if (@_ == 4) {
		($conversion_list) = $_[3];
	}
	
	#Read in the species list.
	my $full_name = 0;
	my $short_name = 1;
	my $acronym_name = 2;
	my %spp1 = %{ &OrthoReD_library::species_lister ($species_list, $short_name, $acronym_name) };
	my %spp2 = %{ &OrthoReD_library::species_lister ($species_list, $acronym_name, $full_name) };
	
	#Selecting the files that will be treated.
	my @in_files;
	my $in_files = `ls $in_folder`;
	@in_files = split(/\n/, $in_files);
	foreach my $in_file (@in_files) {
		$in_file =~ s/\r//sig;
		$in_file =~ s/\n//sig;
		$in_file = "$in_folder/".$in_file;
	}
	
	#Make a hash of hash of array: first-order hash representing one FASTA file, second-order hash representing each species, third order array representing a list of orthologs.
	#generating the first-order hash.
	my %gois;
	foreach my $in_file (@in_files) {
		#identifying the name of the gene.
		my $gene;
		my @gene = split(/\//, $in_file);
		@gene = split(/_/, $gene[1]);
		$gene = $gene[0]."_".$gene[1]."_".$gene[2]."_".$gene[3];
		
		#generating the second-order hash.
		my %spp;
		foreach my $sp (sort keys %spp2) {
			#generating the third-order array.
			my @orthologs;
			$spp{$spp2{$sp}} = \@orthologs;
		}
		$gois{$gene} = \%spp;
	}
	
	#Add the genes into the multi-order hash.
	foreach my $in_file (@in_files) {
		#identifying the name of the gene.
		my $gene;
		my @gene = split(/\//, $in_file);
		@gene = split(/_/, $gene[1]);
		$gene = $gene[0]."_".$gene[1]."_".$gene[2]."_".$gene[3];
		
		#reading in a file with orthologs.
		open (INFILE, "<$in_file") or die "cannot open $in_file.\n";
		while (my $line = <INFILE>) {
			if ($line !~ m/^>/) {
				next;
			}
			
			#identifying the header.
			my $header = $line;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			$header =~ s/>//;
			
			#identifying the species the header belongs to.
			my $sp;
			my @sp = split(/_/, $header);
			$sp = $spp2{$sp[0]};
			
			#adding orthologs to each array.
			my %spp = %{ $gois{$gene} };
			my @orthologs = @{ $spp{$sp} };
			
			push (@orthologs, $header);
			$spp{$sp} = \@orthologs;
			$gois{$gene} = \%spp;
		}
	}
	
	#Convert the multi-order hash into a table format.
	open (TEMP, ">temp.txt") or die "cannot open $out_file.\n";
	my @spp;
	foreach my $sp (sort keys %spp2) {
		push (@spp, $spp2{$sp});
	}
	@spp = sort @spp;
	my $spp = join("\t", @spp);
	$spp = "\t".$spp;
	print TEMP "$spp\n";
	foreach my $goi (sort keys %gois) {
		my $line = "$goi";
		my %spp = %{ $gois{$goi} };
		foreach my $sp (sort keys %spp){
			my @orthologs = @{ $spp{$sp} };
			@orthologs = sort @orthologs;
			my $orthologs =join("," , @orthologs);
			$line = $line."\t".$orthologs;
		}
		print TEMP "$line\n";
	}
	close (TEMP);
	
	#Apply the conversion
	if ($conversion_list =~ m/^NONE$/) {
		system "mv temp.txt $out_file";
	}
	if ($conversion_list !~ m/^NONE$/) {
		#Read in the conversion list.
		my $genbankid = 0;
		my $genomeid = 1;
		my %converter;
		open (CONV, "<$conversion_list") or die "cannot open $conversion_list.\n";
		while (my $conversion = <CONV>) {
			$conversion =~ s/\r//sig;
			$conversion =~ s/\n//sig;
			my @line = split (/\t/, $conversion);
			if ($line[1] !~ m/NA/) {
				$converter{$line[$genomeid]} = $line[$genbankid];
			}
		}
		close (CONV);
		
		#Apply the conversion to the temp.txt.
		open (TEMP, "<temp.txt") or die "cannot open temp.txt.\n";
		open (OUT, ">$out_file") or die "cannot open $out_file.\n";
		while (my $line = <TEMP>) {
			$line =~ s/\r//sig;
			$line =~ s/\n//sig;
			my @line = split(/\t/, $line);
			foreach my $element1 (@line) {
				if ($element1 =~ m/,/) {
					my @element1 = split(/,/, $element1);
					foreach my $element2 (@element1) {
						if (exists $converter{$element2}) {
							$element2 = $converter{$element2};
						}
					}
					$element1 = join(",", @element1);
				}
				if (exists $converter{$element1}) {
					$element1 = $converter{$element1};
				}
			}
			$line = join("\t", @line);
			print OUT "$line\n";
		}
		close (TEMP);
		close (OUT);
		unlink ("temp.txt");
	}
}

#Reads a tabulated BLAST output file and parses out top n hits according to the best e-values.
#When there are splice variants, only the variant with the highest e-value will be kept.
sub top_finder {
	if (@_ != 3) {
		die "USAGE: module top_finder requires four arguments: <\$hit_file> <\$top_file> <\$loci_threshold>.\n";
	}
	
	my ($hit_file, $top_file, $loci_threshold) = @_;
	
	my $query = `head -n 1 $hit_file | cut -f 1`;
	$query =~ s/\r//sig;
	$query =~ s/\n//sig;
	
	#prepareing parameters
	my %seqcounts;
	my %lowest_e;
	my %keepers;
	
	#preparing the file.
	my $top_temp = "TOP_FINDER_01_".$query."\.txt";
	system "sort $hit_file -k3,3g -k5,5nr -k4,4nr > $top_temp";
	
	#selecting subjects.
	open (HIT_FILE, "<$top_temp") or die "cannot open $hit_file.\n";
	open (TOP_FILE, ">$top_file") or die "cannot open $top_file.\n";
	while (my $line = <HIT_FILE>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		
		my @line = split ("\t", $line);
		my $subject = $line[0];
		my $hit = $line[1];
		my $eval = $line[2];
		
		@line = split ("_", $hit);
		my $sp = $line[0];
		my $locus = $line[0]."_".$line[1];
		
		if (exists $keepers{$locus}) {
			next;
		}
		elsif ($eval == 0) {
			if (exists $seqcounts{$sp}) {
				$seqcounts{$sp}++;
			}
			else {
				$seqcounts{$sp} = 1;
			}
			$lowest_e{$sp} = 0;
			$keepers{$locus} = "KEEP";
			print TOP_FILE "$line\n";
		}
		elsif (!exists $seqcounts{$sp}) {
			$seqcounts{$sp} = 1;
			$lowest_e{$sp} = $eval;
			$keepers{$locus} = "KEEP";
			print TOP_FILE "$line\n";
		}
		elsif ($seqcounts{$sp} < $loci_threshold or $lowest_e{$sp} == $eval) {
			$seqcounts{$sp}++;
			$lowest_e{$sp} = $eval;
			$keepers{$locus} = "KEEP";
			print TOP_FILE "$line\n";
		}
		else {
			next;
		}
	}
	close (HIT_FILE);
	close (TOP_FILE);
	unlink ($top_temp);
}

#Sets up the translation from nucleic acid to amino acid.
sub translator1 {
	if (@_ != 0) {
		die "USAGE: module translator1 does not take any arguments.\n";
	}
	my %translator1 = (
		'GCA' => 'A', #Alanine
		'GCB' => 'A', #Alanine
		'GCC' => 'A', #Alanine
		'GCD' => 'A', #Alanine
		'GCG' => 'A', #Alanine
		'GCH' => 'A', #Alanine
		'GCK' => 'A', #Alanine
		'GCM' => 'A', #Alanine
		'GCN' => 'A', #Alanine
		'GCR' => 'A', #Alanine
		'GCS' => 'A', #Alanine
		'GCT' => 'A', #Alanine
		'GCV' => 'A', #Alanine
		'GCW' => 'A', #Alanine
		'GCY' => 'A', #Alanine
		'TGC' => 'C', #Cysteine
		'TGT' => 'C', #Cysteine
		'TGY' => 'C', #Cysteine
		'GAC' => 'D', #Aspartic Acid
		'GAT' => 'D', #Aspartic Acid
		'GAY' => 'D', #Aspartic Acid
		'GAA' => 'E', #Glutamic Acid
		'GAG' => 'E', #Glutamic Acid
		'GAR' => 'E', #Glutamic Acid
		'TTC' => 'F', #Phenylalanine
		'TTT' => 'F', #Phenylalanine
		'TTY' => 'F', #Phenylalanine
		'GGA' => 'G', #Glycine
		'GGB' => 'G', #Glycine
		'GGC' => 'G', #Glycine
		'GGD' => 'G', #Glycine
		'GGG' => 'G', #Glycine
		'GGH' => 'G', #Glycine
		'GGK' => 'G', #Glycine
		'GGM' => 'G', #Glycine
		'GGN' => 'G', #Glycine
		'GGR' => 'G', #Glycine
		'GGS' => 'G', #Glycine
		'GGT' => 'G', #Glycine
		'GGV' => 'G', #Glycine
		'GGW' => 'G', #Glycine
		'GGY' => 'G', #Glycine
		'CAC' => 'H', #Histidine
		'CAT' => 'H', #Histidine
		'CAY' => 'H', #Histidine
		'ATA' => 'I', #Isoleucine
		'ATC' => 'I', #Isoleucine
		'ATH' => 'I', #Isoleucine
		'ATM' => 'I', #Isoleucine
		'ATT' => 'I', #Isoleucine
		'ATW' => 'I', #Isoleucine
		'ATY' => 'I', #Isoleucine
		'AAA' => 'K', #Lysine
		'AAG' => 'K', #Lysine
		'AAR' => 'K', #Lysine
		'CTA' => 'L', #Leucine
		'CTC' => 'L', #Leucine
		'CTD' => 'L', #Leucine
		'CTG' => 'L', #Leucine
		'CTH' => 'L', #Leucine
		'CTK' => 'L', #Leucine
		'CTM' => 'L', #Leucine
		'CTN' => 'L', #Leucine
		'CTR' => 'L', #Leucine
		'CTS' => 'L', #Leucine
		'CTT' => 'L', #Leucine
		'CTV' => 'L', #Leucine
		'CTV' => 'L', #Leucine
		'CTW' => 'L', #Leucine
		'CTY' => 'L', #Leucine
		'TTA' => 'L', #Leucine
		'TTG' => 'L', #Leucine
		'TTR' => 'L', #Leucine
		'YTA' => 'L', #Leucine
		'YTG' => 'L', #Leucine
		'YTR' => 'L', #Leucine
		'ATG' => 'M', #Methionine
		'AAC' => 'N', #Asparagine
		'AAT' => 'N', #Asparagine
		'AAY' => 'N', #Asparagine
		'CCA' => 'P', #Proline
		'CCB' => 'P', #Proline
		'CCC' => 'P', #Proline
		'CCD' => 'P', #Proline
		'CCG' => 'P', #Proline
		'CCH' => 'P', #Proline
		'CCK' => 'P', #Proline
		'CCM' => 'P', #Proline
		'CCN' => 'P', #Proline
		'CCR' => 'P', #Proline
		'CCS' => 'P', #Proline
		'CCT' => 'P', #Proline
		'CCV' => 'P', #Proline
		'CCW' => 'P', #Proline
		'CCY' => 'P', #Proline
		'CAA' => 'Q', #Glutamine
		'CAG' => 'Q', #Glutamine
		'CAR' => 'Q', #Glutamine
		'AGA' => 'R', #Arginine
		'AGG' => 'R', #Arginine
		'AGR' => 'R', #Arginine
		'CGA' => 'R', #Arginine
		'CGB' => 'R', #Arginine
		'CGC' => 'R', #Arginine
		'CGD' => 'R', #Arginine
		'CGG' => 'R', #Arginine
		'CGH' => 'R', #Arginine
		'CGK' => 'R', #Arginine
		'CGM' => 'R', #Arginine
		'CGN' => 'R', #Arginine
		'CGR' => 'R', #Arginine
		'CGS' => 'R', #Arginine
		'CGT' => 'R', #Arginine
		'CGV' => 'R', #Arginine
		'CGW' => 'R', #Arginine
		'CGY' => 'R', #Arginine
		'MGA' => 'R', #Arginine
		'MGG' => 'R', #Arginine
		'MGR' => 'R', #Arginine
		'AGC' => 'S', #Serine
		'AGT' => 'S', #Serine
		'AGY' => 'S', #Serine
		'TCA' => 'S', #Serine
		'TCC' => 'S', #Serine
		'TCD' => 'S', #Serine
		'TCG' => 'S', #Serine
		'TCH' => 'S', #Serine
		'TCK' => 'S', #Serine
		'TCM' => 'S', #Serine
		'TCN' => 'S', #Serine
		'TCR' => 'S', #Serine
		'TCS' => 'S', #Serine
		'TCT' => 'S', #Serine
		'TCV' => 'S', #Serine
		'TCV' => 'S', #Serine
		'TCW' => 'S', #Serine
		'TCY' => 'S', #Serine
		'ACA' => 'T', #Threonine
		'ACC' => 'T', #Threonine
		'ACD' => 'T', #Threonine
		'ACG' => 'T', #Threonine
		'ACH' => 'T', #Threonine
		'ACK' => 'T', #Threonine
		'ACM' => 'T', #Threonine
		'ACN' => 'T', #Threonine
		'ACR' => 'T', #Threonine
		'ACS' => 'T', #Threonine
		'ACT' => 'T', #Threonine
		'ACV' => 'T', #Threonine
		'ACV' => 'T', #Threonine
		'ACW' => 'T', #Threonine
		'ACY' => 'T', #Threonine
		'GTA' => 'V', #Valine
		'GTC' => 'V', #Valine
		'GTD' => 'V', #Valine
		'GTG' => 'V', #Valine
		'GTH' => 'V', #Valine
		'GTK' => 'V', #Valine
		'GTM' => 'V', #Valine
		'GTN' => 'V', #Valine
		'GTR' => 'V', #Valine
		'GTS' => 'V', #Valine
		'GTT' => 'V', #Valine
		'GTV' => 'V', #Valine
		'GTV' => 'V', #Valine
		'GTW' => 'V', #Valine
		'GTY' => 'V', #Valine
		'TGG' => 'W', #Tryptophan
		'TAC' => 'Y', #Tyrosine
		'TAT' => 'Y', #Tyrosine
		'TAY' => 'Y', #Tyrosine
		'TAA' => '*', #Stop
		'TAG' => '*', #Stop
		'TAR' => '*', #Stop
		'TGA' => '*', #Stop
		'TRA' => '*', #Stop
	);
	return (%translator1);
}

#Reads a single-lined DNA FASTA file and generates a translated single-lined AA FASTA file.
#Any codon that cannot be translated to a specific amino acid will be translated as X.
sub translator2 {
	if (@_ != 2) {
		die "USAGE: module translator2 requires three arguments: <\$dna_file> <\$aa_file>.\n";
	}
	
	my ($dna_file, $aa_file) = @_;
	
	my %translator1 = &OrthoReD_library::translator1 ();
	
	open (CDS_DNA, "<$dna_file") or die "cannot open $dna_file.\n";
	open (CDS_AA, ">$aa_file") or die "cannot open $aa_file.\n";
	while (my $line = <CDS_DNA>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		my $headder;
		my $cds;
		my $transcript;
		if ($line =~ m/^>/) {
			$headder = $line;
			$cds = <CDS_DNA>;
			$cds =~ s/\r//sig;
			$cds =~ s/\n//sig;
			$cds =~ s/^CTG/ATG/;
			$cds =~ s/^TTG/ATG/;
			my $cds_length = length ($cds);
			for (my $c = 1; $c <= $cds_length / 3; $c++) {
				my $codon = substr ($cds, $c * 3 - 3, 3);
				if (exists $translator1{$codon}) {
					$transcript = $transcript.$translator1{$codon};
				}
				else {
					$transcript = $transcript."X";
				}
			}
			print CDS_AA "$headder\n";
			print CDS_AA "$transcript\n";
		}
		else {
			next;
		}
	}
	close (CDS_DNA);
	close (CDS_AA);
}

#Given a tree, gene of interest and a threshold, generates a new tree with only taxa connected to the gene of interest with branches under the threshold
sub tree_cutter {
	if (@_ != 4) {
		die "USAGE: module tree_cutter requires four arguments: <\$tree> <\$gene> <\$dist_threshold> <\$out_tre>.\n";
	}

	my ($tree_file, $gene, $dist_threshold, $out_tre) = @_;
	
	#cleaning the tree
	my $tree = `head -n 1 $tree_file`;
	$tree =~ s/\n//sig;
	$tree =~ s/\r//sig;
	$tree =~ s/\;//;
	$tree = $tree."\:0";
	
	my $test_tree = $tree;
	my $cut_test = 1;
	while ($cut_test == 1) {
		#cut one branch too long at a time
		($test_tree, $cut_test) = &long_branch_cutter ($test_tree, $gene, $dist_threshold);
	}
	$tree = $test_tree;
	
	open (OUT, ">$out_tre") or die "cannot open $out_tre.\n";
	print OUT "$tree\;\n";
	close (OUT);
}

#Reads a reduced tree and the conversion list, and generates the expanded tree.
sub tree_expander {
	if (@_ != 3) {
		die "USAGE: module tree_expander requires three arguments: <\$reduced_tre> <\$conv_list> <\$full_tre>.\n";
	}
	
	my ($reduced_tre, $conv_list, $full_tre) = @_;
	
	#Read the conversion file
	my %conversion;
	
	open (CONV, "<$conv_list") or die "cannot open $conv_list.\n";
	while (my $line = <CONV>){
		$line =~ s/\n//sig;
		$line =~ s/\r//sig;
		
		my @list = split (/\t/, $line);
		my $length = scalar (@list);
		
		if ($length == 1) {
			next;
		}
		else {
			my $replace = join ("\:0\,", @list);
			$replace = "\(".$replace."\:0\)";
			
			$conversion{$list[0]} = $replace;
		}
	}
	close (CONV);

	#Read the tre to be converted
	my $tre = `head -n 1 $reduced_tre`;
	$tre =~ s/\n//sig;
	$tre =~ s/\r//sig;
	
	foreach my $convert (sort keys %conversion) {
		$tre =~ s/$convert/$conversion{$convert}/;
	}
	
	#open output file
	open (OUT, ">$full_tre") or die "cannot open $full_tre\n";
	print OUT "$tre\n";
	close (OUT);
}

#Reads a tabulated BLAST output file and parses out unique BLAST hits by keeping the hit with the best raw alignemnt score.
sub uniq_finder {
	if (@_ != 2) {
		die "USAGE: module uniq_finder requires two arguments: <\$all_file> <\$uniq_file>.\n";
	}
	
	my ($hit_file, $uniq_file) = @_;
	
	my $subject;
	my @subject = split (/\//, $hit_file);
	@subject = split (/\./, $subject[-1]);
	@subject = split (/_/, $subject[0]);
	$subject = $subject[0]."_".$subject[1]."_".$subject[2]."_".$subject[3];
	
	my $uniq_temp = "UNIQ_FINDER_01_".$subject."\.txt";
	system "sort -k3,3g $hit_file >$uniq_temp";
	
	my %uniq;
	open (HIT_FILE, "<$uniq_temp") or die "cannot open $uniq_temp.\n";
	open (UNIQ_FILE, ">$uniq_file") or die "cannot open $uniq_file.\n";
	while (my $line = <HIT_FILE>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		
		my @line = split (/\t/, $line);
		my $hit = $line[1];
		
		if (!exists $uniq{$hit}) {
			print UNIQ_FILE "$line\n";
			$uniq{$hit} = "EXISTS";
		}
		else {
			next;
		}
	}
	close (HIT_FILE);
	close (UNIQ_FILE);
	
	unlink ($uniq_temp);
}

#Reads a tabulated BLAST output file and parses out unique BLAST hits by keeping the hit with the best e-value.
sub uniq_pair_finder {
	if (@_ != 2) {
		die "USAGE: module uniq_finder requires two arguments: <\$all_file> <\$uniq_file>.\n";
	}
	
	my ($hit_file, $uniq_file) = @_;
	
	my %uniq;
	open (HIT_FILE, "<$hit_file") or die "cannot open $hit_file.\n";
	while (my $line = <HIT_FILE>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		
		my @line = split (/\t/, $line);
		my $pair = $line[0].$line[1];
		
		if (exists $uniq{$pair}) {
			my @saved_line = split (/\t/, $uniq{$pair});
			
			if ($saved_line[3] < $line[3]) {
				$uniq{$pair} = $line;
			}
			else {
				next;
			}
		}
		if (!exists $uniq{$pair}) {
			$uniq{$pair} = $line;
		}
	}
	close (HIT_FILE);
	
	open (UNIQ_FILE, ">$uniq_file") or die "cannot open $uniq_file.\n";
	foreach my $key (sort keys %uniq) {
		print UNIQ_FILE "$uniq{$key}\n";
	}
	close (UNIQ_FILE);
}

########################################################################
#                                                                      #
#  Following portion is part taken from OrthologID.                    #
#                                                                      #
#  OrthologID is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published   #
#  by the Free Software Foundation, either version 3 of the License,   #
#  or (at your option) any later version.                              #
#                                                                      #
#  OrthologID is distributed in the hope that it will be useful,but    #
#  WITHOUT ANY WARRANTY; without even the implied warranty of          #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    #
#  General Public License for more details.                            #
#                                                                      #
#  You should have received a copy of the GNU General Public License   #
#  along with OrthologID. If not, see <http://www.gnu.org/licenses/>.  #
#                                                                      #
#                                                                      #
#                               Author: Ernest K Lee <elee@amnh.org>   #
#                               Copyright (C) 2006-2011 Ernest K. Lee  #
#                                                                      #
########################################################################

#Return list of leaves of a node belonging to species in arguments
sub getLeaves {
	if (@_ != 2) {
		die "USAGE: module getLeaves requires two arguments: <\$tree_object> <\$spp_list>.\n";
	}
	
	my ($node, $splist) = @_;
	my @sp = @{ $splist };
	my @leaves = ();
	
	foreach my $child (@{ $node->{"child"} }) {
		if (defined $child->{"label"}) {
			my $label = $child->{"label"};
			foreach my $sp (@sp) {
				if ($label =~ /^$sp/) {
					push(@leaves, $label);
					last;
				}
			}
		}
		else {
			push(@leaves, getLeaves($child, $splist));
		}
	}
	return @leaves;
}

#Create new tree node. Returns pointer.
sub newNode{
	if (@_ != 2 and @_ != 3) {
		die "USAGE: module newNode requires two or three arguments: <\$parent_node> <\$node_number> <\$node_label>.\n";
	}
	
	my ($parent, $nodeNum, $label) = @_;
	
	my %node = (
		num => $nodeNum, #Node number in DFS order
		label => $label, #Taxon name
		child => [],
		parent => $parent,
		spMember => {}, #Taxa species of this node (clade)
		eStatus => 1, #Evolution status of node (pure dup/sp[1] or mix[0])
	);
	return \%node;
}

#Return list of ortholog groups as list of list refs
sub orthologGroups {
	if (@_ != 2) {
		die "USAGE: module orthologGroups requires two arguments: <\$tree_object> <\$spp_list>.\n";
	}
	
	my ($treeNode, $orthSp) = @_;
	my @orthSp = @{ $orthSp };
	my @oGroups = ();
	#print "@orthSp\n";
	
	if ($treeNode->{"eStatus"}) {
		my $spCount = 0;
		foreach (@orthSp) {
			$spCount++ if defined($treeNode->{"spMember"}->{$_});
		}
		my @group = ();
		if ($spCount > 1) {
			@group = &getLeaves($treeNode, $orthSp);
		}
		if (@group > 0) {
			@oGroups = (\@group);
		}
	}
	else {
		foreach my $child (@{ $treeNode->{"child"} }) {
			$orthSp = \@orthSp;
			my @groups = &orthologGroups($child, $orthSp);
			if (@groups > 0) {
				push(@oGroups, @groups);
			}
		}
	}
	return @oGroups;
}

#Print tree nodes and their member species (for debugging)
sub printNodes {
	if (@_ != 1) {
		die "USAGE: module printNodes requires one argument: <\$tree_object>.\n";
	}
	
	my $treeObj = shift;
	
	print "node #".$treeObj->{"num"}." [".$treeObj->{"eStatus"}."]";
	while (my($sp, $num) = each(%{ $treeObj->{"spMember"} })) {
		print " $sp: $num ";
	}
	if (defined $treeObj->{"label"}) {
		print $treeObj->{"label"};
	}
	print "\n";
	foreach (@{ $treeObj->{"child"} }) {
		&printNodes($_);
	}
}

# Generate a tree object from a parenthetical tree
sub toTreeObj {
	if (@_ != 3) {
		die "USAGE: module toTreeObj requires three arguments: <\$tre_file> <\$spp_list> <\$overlap_threshold>.\n";
	}
	
	my ($pTree, $spPrefix, $o) = @_;
	
	my $nodeNum = 0;
	my $unknownSp = "Unknown";
	
	#Create root node, remove root parentheses
	my $treeObj = &newNode(undef, $nodeNum++);
	my $currNode;
	
	#Traverse the parenthetical tree
	while ($pTree =~ m/([\(\),])([^\(\),]*)/g) {
		my ($paren, $name) = ($1, $2);
		if ($paren eq "(") {
			if (!defined $currNode) {
				$currNode = $treeObj;
			}
			else {
				my $newNode = &newNode($currNode, $nodeNum++);
				push(@{ $currNode->{"child"} }, $newNode);
				$currNode = $newNode;
			}
		}
		elsif ($paren eq ")") {
			#Calculate evo status of current node
			my %spMem;
			my $spDup = 0; #number of species appear in more than 1 children
			foreach my $ch (@{ $currNode->{"child"} }) {
				$currNode->{"eStatus"} *= $ch->{"eStatus"};
				if ($currNode->{"eStatus"} == 0) {
					last; #nothing to check
				}
				foreach (@$spPrefix) {
					if (defined($ch->{"spMember"}->{$_})) {
						if (defined $spMem{$_}) {
							$spDup++;
						}
						else {
							$spMem{$_} = 1;
						}
					}
				}
			}
			#If more than 1 species present in all children, treat as duplication
			if ($spDup > $o) {
				$currNode->{"eStatus"} = 0;
				
				#NFC-specific
				my @spp = ();
				push (@spp, keys %{ $currNode->{"spMember"} });
				my $spp = \@spp;
				my $test = &order_tester($spp);
				
				if ($test =~ m/YES/) {
					$currNode->{"eStatus"} = 1;
				}
			}
			
			# Pure duplications if only one species
			if (keys(%{ $currNode->{"spMember"} }) == 1) {
				$currNode->{"eStatus"} = 1;
			}
			my $parent = $currNode->{"parent"};
			
			#Done if we are at root
			if (!defined $parent) {
				last;
			}
			
			#Propagate species members to parent
			foreach (@$spPrefix) {
				if (defined $currNode->{"spMember"}->{$_}) {
					$parent->{"spMember"}->{$_} += $currNode->{"spMember"}->{$_};
				}
			}
			
			#Move up
			$currNode = $parent;
		}
		
		if ($name ne "") {
			#Create leaf child (taxon) and add species member
			my $newNode = &newNode($currNode, $nodeNum++, $name);
			push(@{ $currNode->{"child"} }, $newNode);
			my $spMember = $currNode->{"spMember"};
			my $matched = 0;
			foreach (@$spPrefix) {
				if ($name =~ /^$_/) {
					$spMember->{$_}++; #Add taxon species to count
					$newNode->{"spMember"}->{$_} = 1;
					$matched = 1;
					last;
				}
			}
			if (!$matched) { #should not happen
				#Unknown species
				die "Unknown species ($name) encountered in tree!\n";
				$spMember->{$unknownSp}++;
				$newNode->{"spMember"}->{$unknownSp}++;
			}
		}
	}
	return $treeObj;
}

sub order_tester {
	if (@_ != 1) {
		die "USAGE: module order_tester requires one argument: <\$spp>.\n";
	}
	
	my $spp = shift;
	my @spp = @{ $spp };
	my $test = "NO";
	
	my %order = (
		'At' => 'Brassicales',
		'Ct' => 'Rosales',
		'Cf' => 'Fabales',
		'Cs' => 'Cucurbitales',
		'Dg' => 'Cucurbitales',
		'Fv' => 'Rosales',
		'Gm' => 'Fabales',
		'Lj' => 'Fabales',
		'Md' => 'Rosales',
		'Me' => 'Malpighiales',
		'Mt' => 'Fabales',
		'Pv' => 'Fabales',
		'Pt' => 'Malpighiales',
		'Pp' => 'Rosales',
		'Vv' => 'Vitales'
	);
	
	my %present_order;
	foreach my $sp (@spp) {
		if (exists $order{$sp}) {
			if (!exists $present_order{$order{$sp}}) {
				$present_order{$order{$sp}} = "EXISTS";
			}
			else {
				next;
			}
		}
		else {
			next;
		}
	}
	my @norders = sort keys %present_order;
	my $norders = @norders;
	
	if ($norders < 2) {
		$test = "YES";
	}
	
	return $test;
}

1;
__END__
