#/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
use Data::Dumper;

if (@ARGV < 2) {
	warn "
	Description: this program is used to count core, pan, dispensable and specific genes. And make the list file.
	             Input file is the cluster file of software 'cd-hit'.
	Usage: perl stat_list.pl <in cluster file> <in sample ID list> [Options]
	       <in sample ID list>: make sure the first column is sample id.
		   Options:
		   --prefix <str>           set ouput files' prefix. Default 'all'.
		   --outdir <path>          set the output directory path. Default './'.
		   filter options:
		   --min_identity <float>   set the minimal identity of the alignment. Default 0.
		   --min_Rcov <float>       set the minimal percentage of representative sequence's (align length)/(sequence length).
		                            Default 0.
		   --min_Scov <float>       set the minimal percentage of alignment sequence's (align length)/(sequence length).
		                            Default 0.\n";
	exit 1;
}

my ($Cluster, $ID_list) = @ARGV;
foreach ($Cluster, $ID_list) {-s $_ || die "Input file is not exist or empty $_\n";}
my ($Prefix, $Outdir, $Min_identity, $Min_Rcov, $Min_Scov) = ('all', '.', 0, 0, 0);
GetOptions(
	"prefix:s" => \$Prefix, "outdir:s" => \$Outdir, "min_identity:f" => \$Min_identity,
	"min_Rcov:f" => \$Min_Rcov, "min_Scov:f" => \$Min_Scov
);

mkdir $Outdir unless (-d $Outdir);

my (%Sample_ID, @Sample_ID);
open LIST, $ID_list || die;
my $count = 0;
while (<LIST>) {
	chomp;
	next if (/^\s*$|^\s*#/);
	if (/^\s*(\S+)/ && !exists $Sample_ID{$1}) {
		$Sample_ID{$1} = ++$count;
	}
}
close LIST;
@Sample_ID = sort {$Sample_ID{$a} <=> $Sample_ID{$b}} keys %Sample_ID;

open OUT_Pan, ">$Outdir/$Prefix\_Pan.matrix" || die;
open OUT_Core, ">$Outdir/$Prefix\_Core.matrix" || die;
open OUT_Disp, ">$Outdir/$Prefix\_Dispensable.matrix" || die;
open OUT_Pan_v, ">$Outdir/$Prefix\_Pan.venn" || die;
open OUT_Cluster, ">$Outdir/$Prefix\_GeneCluster.list" || die;
print OUT_Pan join ("\t", "", @Sample_ID), "\n";
print OUT_Core join ("\t", "", @Sample_ID), "\n";
print OUT_Disp join ("\t", "", @Sample_ID), "\n";
print OUT_Pan_v join ("\t", "", @Sample_ID), "\n";
print OUT_Cluster "#Pan_Genes\tGene_number\tGene_list\n";
my %Specific;
my ($Pan_number, $Core_number, $Disp_number, $all_number);
open IN, "$Cluster" || die;
$/ = "\n\>";
while (<IN>) {
	chomp;
	my @Pan_Gene; 
	my (%hash_identity, %hash_gene);
#/\n\d+\s+(\d+)aa,\s+\>((\S+)\_([^_\s]+))\.{3}\s+\*/s || next;
    /\n\d+\s+(\d+)aa,\s+\>((\S+)__([^\s]+))\.{3}\s+\*/s || next;
	push @Pan_Gene, "$3($4)";
	push @{$hash_identity{$4}}, 100;
	push @{$hash_gene{$4}}, $3;
	$all_number += 1;
	my ($rep_length, $rep_spe_id) = ($1, $4);
#while (s/\n\d+\s+(\d+)aa,\s+\>((\S+)\_([^_\s]+))\.{3}\s+at\s+(\d+):(\d+):(\d+):(\d+)\/([\d\.]+)\%//s) {
    while (s/\n\d+\s+(\d+)aa,\s+\>((\S+)__([^\s]+))\.{3}\s+at\s+(\d+):(\d+):(\d+):(\d+)\/([\d\.]+)\%//s) {  
    	$all_number += 1;
		($9 < $Min_identity*100 || ($8-$7+1)/$rep_length < $Min_Rcov || ($6-$5+1)/$1 < $Min_Scov) && next;
		push @Pan_Gene, "$3($4)";
		push @{$hash_identity{$4}}, $9;
		push @{$hash_gene{$4}}, $3;
	}
	print OUT_Cluster "$Pan_Gene[0]\t[$#Pan_Gene]\t", join (" ", sort @Pan_Gene[1..$#Pan_Gene]), "\n";
	my $spe_number = 0;
	my $out_str = "$Pan_Gene[0]";
	my $out_venn = "$Pan_Gene[0]";
	$Pan_number ++;
	foreach (@Sample_ID) {
		if (exists $hash_identity{$_}) {
			$spe_number ++;
			$out_str .= "\t" . (sort {$b <=> $a} @{$hash_identity{$_}})[0];
			$out_venn .= "\t$Pan_number";
		}
		else {
			$out_str .= "\t0";
			$out_venn .= "\t-";
		}
	}
	print OUT_Pan "$out_str\n";
	print OUT_Pan_v "$out_venn\n";
	if ($spe_number == @Sample_ID) { print OUT_Core "$out_str\n"; $Core_number ++;}
	else { print OUT_Disp "$out_str\n"; $Disp_number ++;}
	if ((my ($id) = keys %hash_identity) == 1) {
		push @{$Specific{$id}}, @{$hash_gene{$id}};
	}
}
$/ = "\n";
close IN;
close OUT_Cluster;
close OUT_Pan_v;
close OUT_Disp;
close OUT_Core;
close OUT_Pan;

open OUT_Spec, ">$Outdir/$Prefix\_Specific.list" || die;
print OUT_Spec "#Specific_name\tGene_number\tGene_list\n";
my $Specific_number;
foreach (@Sample_ID) {
	if (exists $Specific{$_}) {
		my $gene_number = @{$Specific{$_}};
		print OUT_Spec "$_\t$gene_number\t", join (" ", sort @{$Specific{$_}}), "\n";
		$Specific_number .= "\t$_\t$gene_number\n";
	}
	else { print OUT_Spec "$_\t0\n"; $Specific_number .= "\t$_\t0\n";}
}
close OUT_Spec;

#------------- output stat info
open OUT_STAT, ">$Outdir/$Prefix\_Stat.xls" || die;
print OUT_STAT
	"All Gene(#):\t$all_number\n",
	"Pan Gene(#):\t$Pan_number\n",
	"Core Gene(#):\t$Core_number\n",
	"Dispensable Gene(#):\t$Disp_number\n",
	"\nStrain Specific Gene:\n",
	"\tSpecies ID\tNumber(#):\n",
	"$Specific_number\n";
close OUT_STAT;
__END__

>Cluster 8
0       2297aa, >FGM003568_F... *
>Cluster 9
0       2084aa, >FGM001247_F... at 1:2084:1:2081/79.17%
1       2241aa, >NSU_1236_N.pent-US6-1... *
2       2082aa, >PP1Y_AT6484_N.PP1Y... at 1:2082:3:2081/79.78%
>Cluster 473
0       853aa, >FGM002737_F... at 1:853:26:878/83.35%
1       874aa, >FGM000428_F... at 1:874:5:878/98.51%
2       878aa, >NSU_3777_N.pent-US6-1... *
3       874aa, >NSU_0592_N.pent-US6-1... at 1:874:5:878/97.94%
4       845aa, >PP1Y_Mpl6117_N.PP1Y... at 1:845:34:878/83.55%
>Cluster 581
0       832aa, >FGM000905_F... at 6:832:1:841/75.12%
1       832aa, >FGM005146_F... at 1:832:14:841/82.21%
2       841aa, >FGM000874_F... *

