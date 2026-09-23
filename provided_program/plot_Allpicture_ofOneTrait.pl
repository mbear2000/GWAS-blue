use strict;
use warnings;

if(@ARGV!=3){
	print "\nError:  perl plot_Allpicture_ofOneTrait.pl EMMAx.Result/hIBS/2nd3rd_20GrainWeight2012_Ave_chr01.hIBS 2nd3rd_20GrainWeight2012_Ave 2nd3rd_20GrainWeight2012_Ave\n";
	print "!!!check reference file: /data9/home/yzhao/reference/IRGSP-1.0_genome.fasta.fai		\n\n";
	exit;
}

my $pvalueFile=$ARGV[0];
my $outFile=$ARGV[1];
my $printText=$ARGV[2];

###	perl program/plot_manhattan.pl EMMAx.Result/hIBS/2nd3rd_20GrainWeight2012_Ave_chr01.hIBS 2nd3rd_20GrainWeight2012_Ave 2nd3rd_20GrainWeight2012_Ave
my $text1="/usr/bin/perl program/plot_manhattan.pl $pvalueFile $outFile $printText";
print "$text1\n";
system $text1;

###	perl program/plot_manhattan_eachChr.pl EMMAx.Result/hIBS/2nd3rd_20GrainWeight2012_Ave_chr01.hIBS 2nd3rd_20GrainWeight2012_Ave 2nd3rd_20GrainWeight2012_Ave
my $text2="/usr/bin/perl program/plot_manhattan_eachChr.pl $pvalueFile $outFile $printText";
print "$text2\n";
system $text2;

###	perl program/plot_QQplot2023.pl EMMAx.Result/hIBS/2nd3rd_20GrainWeight2012_Ave_chr01.hIBS 2nd3rd_20GrainWeight2012_Ave 2nd3rd_20GrainWeight2012_Ave
my $text3="/usr/bin/perl program/plot_QQplot2023.pl $pvalueFile $outFile $printText";
print "$text3\n";
system $text3;

