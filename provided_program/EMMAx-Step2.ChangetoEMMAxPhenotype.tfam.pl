use strict;
use warnings;

#change phenotype file to tfam file format which is used in the EMMAx program
#perl EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl phenotype.1.O_rufipogon.446lines.csv O_rufipogon.446lines
my(@value,@sample,$m,@pheno);
my $title=`head -1 $ARGV[0]`;
chomp($title);
my @phenoTab=split /,/,$title;
my $col=@phenoTab;
open INPUT,"<$ARGV[0]" or die "Can't open the INPUT file (Original Phenotype)...\n";
my $n=0;
while(<INPUT>){
	if(!($_=~/^Accession/)){
		$m=1;
		while($m<=$col){
			$_=~s/,,/,NA,/;
			$m++;
		}
		$_=~s/,\n/,NA\n/;
		chomp;
		@value=split /,/,$_;
		$sample[$n]=$value[0];
		$m=1;
		while($m<$col){
			$pheno[$n][$m]=$value[$m];
			$m++;
		}
		$n++;
	}
}
close INPUT;
my $row=$n-1;


$m=1;
while($m<$col){
	open OUTPUT,">EMMAx.Data/$ARGV[1]_$phenoTab[$m].tfam";
	$n=0;
	while($n<=$row){
		print OUTPUT "$sample[$n]\t$sample[$n]\t0\t0\t3\t$pheno[$n][$m]\n";
		$n++;
	}
	close OUTPUT;
	$n=1;
	while($n<=9){
		system "cp tped/$ARGV[1]_chr0$n.tped EMMAx.Data/$ARGV[1]_$phenoTab[$m]_chr0$n.tped";
		system "cp EMMAx.Data/$ARGV[1]_$phenoTab[$m].tfam EMMAx.Data/$ARGV[1]_$phenoTab[$m]_chr0$n.tfam";
		$n++;
	}
	$n=10;
	while($n<=12){
		system "cp tped/$ARGV[1]_chr$n.tped EMMAx.Data/$ARGV[1]_$phenoTab[$m]_chr$n.tped";
		system "cp EMMAx.Data/$ARGV[1]_$phenoTab[$m].tfam EMMAx.Data/$ARGV[1]_$phenoTab[$m]_chr$n.tfam";
		$n++;
	}
	$m++;
}

$m=1;
while($m<$col){
	open OUTPUT,">EMMAx.Data/$ARGV[1]_$phenoTab[$m].pheno";
	$n=0;
	while($n<=$row){
		print OUTPUT "$sample[$n]\t$sample[$n]\t$pheno[$n][$m]\n";
		$n++;
	}
	close OUTPUT;
	$m++;
}
print "finish...\n";
