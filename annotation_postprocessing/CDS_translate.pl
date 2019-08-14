#!/usr/bin/perl

=head1 NAME
CDS_translate.pl

=head1 SYNOPSIS
CDS_translate.pl -i [CDS fasta file] > protein.fasta
script to translate CDS from reading frame offset 0 (from first base) using bioperl based codon usage.
MAKER generated annotations are based on bioperl codon usage.

=head1 COMMAND-LINE OPTIONS
 -i  CDS fasta (required)
 -h  Help

=cut

use strict;
use warnings;
use Bio::SeqIO;
use File::Slurp;
use Getopt::Std;

our ( $opt_i, $opt_h );
getopts('i:h');
if ($opt_h) {
	help();
	exit;
}

my $cds_file = $opt_i;
my $input_cds = read_file($cds_file)
	or die "Could not open input CDS fasta file: $cds_file\n";

my $fasta_file = Bio::SeqIO->new(
	-file   => $cds_file,
	-format => "fasta",
	);

while ( my $cds = $sequences->next_seq ){
	my $protein = $cds->translate(
		-codontable_id => 1, 
		-frame         => 0, 
		);
	print ">", $cds->display_id, "\n";
	print $protein->seq, "\n";
}


sub help {
	print STDERR <<EOF;
	$0:
		Description:
			script to translate CDS from reading frame offset 0 (from first base) using bioperl based codon usage. MAKER generated annotations are based on bioperl codon usage.

		usage:
			CDS_translate.pl -i [CDS fasta file] > protein.fasta

		Flags:
			-i CDS fasta file for tranlation
			-h Help

EOF
	exit(1);
}

__END__