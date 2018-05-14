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

package header_cleaner;

#####MODULES#####
#Reads a header and cleans the header.
sub header_cleaner {
	if (@_ != 2) {
		die "USAGE: module header_cleaner requires two argument: <\$species_list> <\$header>.\n";
	}
	
	my ($species_list, $header) = @_;
	
	my $full_name = 0;
	my $short_name = 1;
	my $acronym_name = 2;
	
	my %spp1 = %{ &OrthoReD_library::species_lister ($species_list, $short_name, $acronym_name) };
	my %spp2 = %{ &OrthoReD_library::species_lister ($species_list, $acronym_name, $full_name) };
	
	#removing ">" at the begining.
	if ($header =~ m/^>/) {
		$header =~ s/^>//;
	}
	
	my $cleaned_header;
	
	#####PLANT DATASET#####
	
	#
	#Arabidopsis thaliana TAIR10
	#
	if ($header =~ m/^AT/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">At_".$field23[0]."_".$field23[1]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Chamaecrista fasciculata v1.1
	#
	if ($header =~/^Cf/ and $header !~ m/GI/) {
		$header =~ s/^Cf//;
		$cleaned_header = ">Cf_".$header."_na_na";
		return ($cleaned_header);
		last;
	}
	
	#
	#Cucumis sativus v1.0
	#
	if ($header =~ m/^Cucsa/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Cs_".$field23[0].$field23[1]."_".$field23[2]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Datisca glomerata transcriptome NR
	#
	if ($header =~ m/^Dg_DgTrNR/) {
		$cleaned_header = ">".$header;
		return ($cleaned_header);
		last;
	}
	
	#
	#Fragaria vesca v1.1
	#
	if ($header =~ m/^mrna/ and $header =~ m/hybrid/ and $header =~ m/pacid/) {
		$header =~ s/://;
		my @header = split (/ /, $header);
		
		my $field23 = $header[0];
		my @field23 = split (/-v1/, $header[0]);
		$field23 = $field23[0];
		@field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Fv_".$field23[0]."_".$field23[1]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Glycine max v1
	#
	if ($header =~ m/^Glyma/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Gm_".$field23[0].$field23[1]."_".$field23[2]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Lotus japonicus v3.0
	#
	if ($header =~ m/^Ljchloro/ or $header =~ m/^Ljmito/ or $header =~ m/g3v/) {
		my @header = split (/ /, $header);
		my @field23 = split (/\./, $header[0]);
		$cleaned_header = ">Lj_".$field23[0]."_".$field23[1]."_na";
		return ($cleaned_header);
		last;
	}
	
	#
	#Mallus domestica
	#
	if ($header =~ m/^MDP/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field2 = $header[0];
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Md_".$field2."_na_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Manihot esculenta
	#
	if ($header =~ m/^Manes/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Me_".$field23[0].$field23[1]."_".$field23[2]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Medicago truncatula
	#
	if ($header =~ m/^Medtr/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Mt_".$field23[0]."_".$field23[1]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Phaseolus vulgaris
	#
	if ($header =~ m/^Phvul/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Pv_".$field23[0].$field23[1]."_".$field23[2]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Populus trichocarpa
	#
	if ($header =~ m/^Potri/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Pt_".$field23[0].$field23[1]."_".$field23[2]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Prunus persica
	#
	if ($header =~ m/^Prupe/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field23 = $header[0];
		my @field23 = split (/\./, $field23);
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Pp_".$field23[0].$field23[1]."_".$field23[2]."_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	#
	#Vitis vinifera
	#
	if ($header =~ m/^GSVIVT/ and $header =~ m/pacid/) {
		my @header = split (/ /, $header);
		my $field2 = $header[0];
		my $field4 = $header[1];
		$field4 =~ s/pacid=//;
		$cleaned_header = ">Vv_".$field2."_na_PACID".$field4;
		return ($cleaned_header);
		last;
	}
	
	##########
	
=cut
	
	#####FLY DATASET#####
	
	#
	#Diptera
	#
	if ($header =~ m/^D/ or $header =~ m/^L/) {
		my @elements = split (/ /, $header);
		my $sp = $elements[0];
		my $gene_id = $elements[1];
		
		$cleaned_header = ">".$sp."_".$gene_id."_na_na";
		return ($cleaned_header);
		last;
	}
	##########
	
=cut
	
	$header = ">".$header;
	die "ERROR: $header could not be cleaned.\n";
}

1;
__END__
