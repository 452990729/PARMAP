#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV != 3) {
	warn "
		Usage: perl specific_gene.pl <in specific list> <in gene file list> <outdir>
		<in specific list> format: species_id gene_number gene_id_list
		<in gene file list> format: species_id gene_file_path\n\n";
	exit;
}

my ($list, $path, $outdir) = @ARGV;
mkdir $outdir unless (-d $outdir);

my (%hash);
open IN, $list || die;
while (<IN>) {
	chomp;
	next if (/^\s*$|^\s*#/);
	s/^\s+|\s+$//g;
	my @w=split;
	($w[1] =~ /\d+/) || next;
	($w[1] > 0) && ($hash{$w[0]} = join " ", ("", @w[2..$#w]), "");
}
close IN;

open IN, $path || die;
while (<IN>) {
	chomp;
	next if (/^\s*$|^\s*#|^\s*\t/);
	s/^\s+|\s+$//g;
	(/^(\S+)\s+(\S+)/) || next;
	my ($species_id, $gene_file) = ($1, $2);
	next if (!exists $hash{$species_id});
	if (!(-e $gene_file && -s $gene_file)) {
		warn "gene file of sample \"$species_id\" is error: $gene_file!\nIgnore it.\n";
		next;
	}
	open OUT, ">$outdir/$species_id\_Specific.fa" || die;
	open GENE, $gene_file || die;
	$/ = "\>"; <GENE>;
	while (<GENE>) {
		chomp;
		/^(\S+)/;
        my $geneID = $1;
        $geneID =~ s/\|//g;
        $hash{$species_id} =~ s/\|//g;
		if ($hash{$species_id} =~ / $geneID /) { print OUT "\>$_";}
	}
	$/ = "\n";
	close GENE;
	close OUT;
}
close IN;







