#!/usr/bin/perl
#	perl plot_QQplot2023.pl <$input> <$input> <$output>
use strict;
use warnings;
use GD;
#Y-observed	X-expected	Ref-maxY
if(@ARGV!=3){
	print "\nError:  perl plot_QQplot2023.pl 2nd3rd_20GrainWeight2012_Ave_chr01.hIBS 2nd3rd_20GrainWeight2012_Ave_hIBS 2nd3rd_20GrainWeight2012_Ave\n\n";
	exit;
}

unless(-e "QQplot"){
	my $dir="mkdir QQplot";
	system $dir;
}

my($Yaxis,$Xaxis,$extend,$radius,$lineWidth,$scale,$fontSize);
################################################
###²ÎÊý
$Yaxis=600;
$Xaxis=600;
$radius=3;
$extend=80;
$lineWidth=3;
$fontSize=25;
################################################

#read in
my $chr=1;
my @pvalue;
my $n=0;
while($chr<=12){
	my $file_input=$ARGV[0];
	if($chr<10){
		$file_input=~s/chr01/chr0$chr/;
	}else{
		$file_input=~s/chr01/chr$chr/;
	}
	open INPUT,"<$file_input" or die "Can't open the INPUT file: $chr\n";
	while(<INPUT>)	{
		chomp;
		if(!($_=~/^CHROM/)){
			my @value=split /\s+/,$_;
			$pvalue[$n][0]=$chr;
			$pvalue[$n][1]=$value[1];
			$pvalue[$n][2]=$value[2];
			$n++;
		}
	}
	close INPUT;
	$chr++;
}
#output
my @sort_pvalue=sort{$b->[2] <=> $a->[2]}@pvalue;
my $number=@sort_pvalue;
my @observed;
my @expected;

#print "$number\n$sort_pvalue[0][0]\t$sort_pvalue[0][1]\t$sort_pvalue[0][2]\n";
#print "$number\n$sort_pvalue[-1][0]\t$sort_pvalue[-1][1]\t$sort_pvalue[-1][2]\n";

$n=0;
open QQoutput,">QQplot/$ARGV[1].value" or die "Can't open the QQoutput file: \n";
while($n<$number){
	$observed[$n]=$sort_pvalue[$n][2];
	$expected[$n]=-1*(log(($n+1)/$number)/log(10)) ;
	print QQoutput "scaffold$sort_pvalue[$n][0]\t$sort_pvalue[$n][1]\t$observed[$n]\t$expected[$n]\n";
	$n++
}
close QQoutput;

#plot
my $timesY=$Yaxis/(int($observed[0])+1);
my $timesX=$Xaxis/(int($expected[0])+1);
#my @change_observed=@observed*100;
#my @change_expected=@expected*100;

open IMAGE,">QQplot/$ARGV[1]_QQplot.png" or die "Can't open the IMAGE file: \n";
my $image=GD::Image->new($Xaxis+$extend*2,$Yaxis+$extend*2);
my $white=$image->colorAllocate(255,255,255);
my $black=$image->colorAllocate(0,0,0);
my $red=$image->colorAllocate(255,0,0);
my $orange=$image->colorAllocate(255,127,0);
my $green=$image->colorAllocate(85,107,47);
my $grey=$image->colorAllocate(192,192,192);
my $blue=$image->colorAllocate(51,102,255);

###four points of rectangle:
#($extend,$extend)							($Xaxis+$extend,$extend)
#($extend,$Yaxis+50)		($Xaxis+$extend,$Yaxis+$extend)

my $width=$Xaxis+$extend;
my $height=$Yaxis+$extend;
#($extend,$extend)				($width,$extend)
#($extend,$height)		($width,$height)

#print "$Xaxis\t$Yaxis\n$width\t$height\n";

#reference line
$n=0;
while($n*$timesX<=$width-$extend)
{
	my $x=$extend+$n*$timesX;
	my $y=$height-$n*$timesY;
	$image->arc($x,$y,$radius,$radius,1,360,$red);
	$image->filledArc($x,$y,$radius,$radius,1,360,$red);
	$n+=0.0001;
}
#maxY

#observed line
$n=0;
while($n<$number)
{
	my $x=$extend+$expected[$n]*$timesX;
	my $y=$height-$observed[$n]*$timesY;
	$image->arc($x,$y,$radius,$radius,1,360,$blue);
	$image->filledArc($x,$y,$radius,$radius,1,360,$blue);
	$n++;
}

##x-axis
$image->rectangle($extend-$lineWidth/2,$height-$lineWidth/2,$width+$lineWidth/2,$height+$lineWidth/2,$black);
$image->filledRectangle($extend-$lineWidth/2,$height-$lineWidth/2,$width+$lineWidth/2,$height+$lineWidth/2,$black);
my $label=0;
while($label*$timesX<=$width-$extend){
	$image->rectangle($extend-$lineWidth/2+$label*$timesX,$height,$extend+$lineWidth/2+$label*$timesX,$height+$lineWidth*3,$black);
	$image->filledRectangle($extend-$lineWidth/2+$label*$timesX,$height,$extend+$lineWidth/2+$label*$timesX,$height+$lineWidth*3,$black);
	$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend+$label*$timesX+3*$lineWidth/2-$fontSize/2,$height+$lineWidth*4+$fontSize,"$label");
	$label++;
}

##y-axis
$image->rectangle($extend-$lineWidth/2,$extend-20,$extend+$lineWidth/2,$height+$lineWidth/2,$black);
$image->filledRectangle($extend-$lineWidth/2,$extend-20,$extend+$lineWidth/2,$height+$lineWidth/2,$black);

$label=0;
my $add;
if(int($observed[0])+1<7){
	$add=1;
}
elsif(int($observed[0])+1>=7 && int($observed[0]+1)<12){
	$add=2;
}
elsif(int($observed[0])+1>=12 && int($observed[0]+1)<30){
	$add=5;
}
elsif(int($observed[0])+1>=30 && int($observed[0]+1)<60){
	$add=10;
}
elsif(int($observed[0])+1>=60 && int($observed[0]+1)<100){
	$add=20;
}
else{
	$add=30;
}
my $extPos;
while($label*$timesY<=$height-$extend)
{
	$image->rectangle($extend-$lineWidth*3,$height-$label*$timesY-$lineWidth/2,$extend,$height-$label*$timesY+$lineWidth/2,$black);
	$image->filledRectangle($extend-$lineWidth*3,$height-$label*$timesY-$lineWidth/2,$extend,$height-$label*$timesY+$lineWidth/2,$black);
	if($label<10){
		$extPos=5;
	}
	elsif($label>=10 && $label<100){
		$extPos=22;
	}
	else{
		$extPos=39;
	}
	$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend-$lineWidth*3-$fontSize/2-$extPos,$height-$label*$timesY+$fontSize/2,"$label");
	$label+=$add;
}
$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",15,0,5,30,"$ARGV[2]");
binmode(IMAGE);
print IMAGE $image->png;
close IMAGE or die "Cant't close IMAGE: $!";
print "QQplot is finished...\n";
