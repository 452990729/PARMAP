#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV != 3) {
	warn "
	Usage: perl core_pan_gene.pl <in gene list> <in gene fa> <out name>\n";
	exit;
}

my ($list, $gene, $out) = @ARGV;
my $id_list = " ";
open IN, $list || die;
while (<IN>) {
	chomp;
	next if (/^\s*$|^\s*#|^\s*\t/);
	$id_list .= "$1\__$2 " if (/^\s*([^\(]+)\((\S+)\)/);
}
close IN;
open OUT, ">$out" || die;
open GENE, $gene || die;
$/ = "\>"; <GENE>;
$/ = "\n";
while (my $head = <GENE>) {
	chomp $head;
	my $id = $1 if ($head =~ /^(\S+)/);
	$/ = "\>";
	chomp (my $seq = <GENE>);
	$/ = "\n";
    $id =~ s/\|//g;
    $id_list =~ s/\|//g;
    if($id_list =~ /$id/){
		print OUT ">$head\n$seq";
	}
}
close GENE;
close OUT;
