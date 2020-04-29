#!/usr/bin/perl -w
use strict;
(@ARGV == 3) || die "perl $0 <pep.list> <out> <min_len>\n";

my ($in, $out, $min_len) = @ARGV;
$min_len ||= 0;
open OUT, ">$out" || die;
open IN, $in || die;
while (<IN>) {
	chomp;
	next if (/^\s*$|^\s*#/);
	s/^\s*|\s*$//g;
	my ($Sid, $path) = ($1, $2) if (/^(\S+)\s+(\S+)/);
	($Sid && -e $path) || ((warn "sample '$Sid' error!\n") && next);
	open SEQ, $path || die;
	$/ = '>'; <SEQ>;
	while (<SEQ>) {
		chomp;
		my @l = split /\n/, $_, 2;
        my $gene_id = (split /\s+/, $l[0])[0];
		$l[1] =~ s/\n+|\s+//g;
		length($l[1]) < $min_len && next;
		print OUT ">$gene_id\__$Sid\n$l[1]\n";
	}
	close SEQ;
	$/ = "\n";
}
close IN;
close OUT;
