use strict;
use warnings;
use GD;

### 阈值=-log((1/snpNumber)*0.01)
##perl plot_manhattan.pl 2nd3rd_20GrainWeight2012_Ave_chr01.hIBS chromosome/2nd3rd_20GrainWeight2012_Ave_hIBS 2nd3rd_20GrainWeight2012_Ave
## input Xaxis Yaxis radius
if(@ARGV!=3){
	print "\nError:  perl plot_manhattan.pl 2nd3rd_20GrainWeight2012_Ave_chr01.hIBS 2nd3rd_20GrainWeight2012_Ave_hIBS 2nd3rd_20GrainWeight2012_Ave\n";
	print "!!!check reference file: /data9/home/yzhao/reference/IRGSP-1.0_genome.fasta.fai		\n\n";
	exit;
}

unless(-e "chromosome"){
	my $dir="mkdir chromosome";
	system $dir;
}

my $chr;
my (@length,%geno_pheno,$position,$snp_num,@maxX,$lower,@value,$Xaxis,$Yaxis,$radius,$extend,$fontSize,$scale,$threshold,$lineWidth);
my $reference="/data9/home/yzhao/reference/IRGSP-1.0_genome.fasta.fai";
open REF,"<$reference" or die "Can't open the reference file: $reference\n\n";
$maxX[0]=0;
while(<REF>){
	chomp;
	@value=split /\t/,$_;
	$chr=$value[0];
	$chr=~s/chr//;
	$chr=$chr+1-1;
	$length[$chr]=$value[1];
	$maxX[$chr]=$value[1]/1000000;
}
close REF;
#######################################
##修改参数
$Xaxis=2000;
$Yaxis=600;
$radius=10;
$extend=80;
$fontSize=25;
#$scale=5;
$threshold=5;
$lineWidth=3;
#######################################

$chr=1;
my $maxY=0;
while($chr<=12){
	my $fileName=$ARGV[0];
	if($chr<10){
		$fileName=~s/chr01/chr0$chr/;
	}else{
		$fileName=~s/chr01/chr$chr/;
	}
	open INPUT, "$fileName" or die "Can't open the INPUT file: $chr";
	#$snp_num=0;
	while(<INPUT>){
		chomp;
		if(!($_=~/^CHROM/)){
			@value=split /\s+/,$_;
			$position=$value[1]/1000000;
			$geno_pheno{$chr}{$value[1]}{0}=$position;
			$geno_pheno{$chr}{$value[1]}{1}=$value[2];
			if($value[2]>$maxY)
			{
				$maxY=$value[2];
			}
			#$snp_num++;
		}
	}
	close INPUT;
	#print "$position\n";
	$chr++;
}

open IMAGE,">chromosome/$ARGV[1].png" or die "Can't open the IMAGE file: $!";
$chr=1;
my $snpLength=0;
while($chr<=12)
{
	$snpLength=$snpLength+$maxX[$chr];
	$chr++;
}
$maxY=int($maxY)+1;
if($maxY<$threshold){
	$maxY=$threshold;
}
my $timesX=$Xaxis/$snpLength;
my $timesY=$Yaxis/$maxY;

if($maxY<=40){
	$scale=5;
}
elsif($maxY>40 && $maxY<=80){
	$scale=20;
}
elsif($maxY>80 && $maxY<=120){
	$scale=30;
}
else{
	$scale=50;
}

my $image=GD::Image->new($Xaxis+$extend*2,$Yaxis+$extend*2);
my $white=$image->colorAllocate(255,255,255);
my $black=$image->colorAllocate(0,0,0);
my $red=$image->colorAllocate(255,0,0);
my $yellow=$image->colorAllocate(255,255,0);
my $orange=$image->colorAllocate(255,127,0);
my $green=$image->colorAllocate(85,107,47);
my $grey=$image->colorAllocate(192,192,192);
my $blue=$image->colorAllocate(51,102,255);

my $color1=$image->colorAllocate(215,0,26);
my $color2=$image->colorAllocate(135,185,52);
my $color3=$image->colorAllocate(252,201,30);
my $color4=$image->colorAllocate(74,33,151);
my $color5=$image->colorAllocate(91,92,24);
my $color6=$image->colorAllocate(214,4,104);
my $color7=$image->colorAllocate(229,113,12);
my $color8=$image->colorAllocate(0,150,211);
my $color9=$image->colorAllocate(215,154,185);
my $color10=$image->colorAllocate(143,205,184);
my $color11=$image->colorAllocate(133,27,41);
my $color12=$image->colorAllocate(12,116,93);

###four points of rectangle:
#($extend,$extend)						($Xaxis+$extend,$extend)
#($extend,$yaxis+$extend)			($Xaxis+$extend,$yaxis+$extend)

my $width=$Xaxis+$extend;
my $height=$Yaxis+$extend;
#($extend,$extend)				($width,$extend)
#($extend,$height)				($width,$height)

##split bar
$chr=1;
$lower=0;
while($chr<12)
{
	$lower=$lower+$maxX[$chr-1];
	$image->rectangle($extend-2+($lower+$maxX[$chr])*$timesX,$extend-20,$extend+($lower+$maxX[$chr])*$timesX,$height,$grey);
	$image->filledRectangle($extend-2+($lower+$maxX[$chr])*$timesX,$extend-20,$extend+($lower+$maxX[$chr])*$timesX,$height,$grey);
	$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,($lower+$maxX[$chr]/2)*$timesX+$extend-$fontSize/2,$height+$fontSize+$lineWidth,"$chr");
	$chr++;
}
$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,($lower+$maxX[8]+$maxX[9]/2)*$timesX+$extend-$fontSize/2,$height+$fontSize+$lineWidth,"12");

##peak
my $num;
my %color;
$color{1}=$color1;
$color{2}=$color2;
$color{3}=$color3;
$color{4}=$color4;
$color{5}=$color5;
$color{6}=$color6;
$color{7}=$color7;
$color{8}=$color8;
$color{9}=$color9;
$color{10}=$color10;
$color{11}=$color11;
$color{12}=$color12;

$chr=1;
$lower=0;
my $the_color;
while($chr<=12){
	$the_color=$color{$chr};
	$lower=$lower+$maxX[$chr-1];
	foreach my $key(sort{$a<=>$b} keys %{$geno_pheno{$chr}}){
		if(exists $geno_pheno{$chr}{$key}{0}){
			my $x=$extend+($geno_pheno{$chr}{$key}{0}+$lower)*$timesX;
			my $y=$height-$geno_pheno{$chr}{$key}{1}*$timesY;
			#print "$x	$y\n";
			$image->arc($x,$y,$radius/2,$radius/2,0,360,$the_color);
			$image->filledArc($x,$y,$radius/2,$radius/2,0,360,$the_color);
		}
	}
	$chr++;
}

##x-axis
$image->rectangle($extend-10-$lineWidth,$height,$width,$height+$lineWidth,$black);
$image->filledRectangle($extend-10-$lineWidth,$height,$width,$height+$lineWidth,$black);
##y-axis
$image->rectangle($extend-$lineWidth,$extend-20,$extend,$height+$lineWidth,$black);
$image->filledRectangle($extend-$lineWidth,$extend-20,$extend,$height+$lineWidth,$black);

my $label=0;
while($label*$timesY<$height-$extend)
{
	$image->rectangle($extend-10,$height-$label*$timesY,$extend,$height-$label*$timesY+$lineWidth,$black);
	$image->filledRectangle($extend-10,$height-$label*$timesY,$extend,$height-$label*$timesY+$lineWidth,$black);
	if($label<10){
		$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend-30,$height-$label*$timesY+$fontSize/2,"$label");
	}
	elsif($label>=10 && $label<100){
		$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend-45,$height-$label*$timesY+$fontSize/2,"$label");
	}
	else{
		$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend-60,$height-$label*$timesY+$fontSize/2,"$label");
	}
	$label+=$scale;
}

$label=$extend;
while($label<$width){
	$image->rectangle($label,$height-$threshold*$timesY-$lineWidth/2,$label+5*$lineWidth/2,$height-$threshold*$timesY+$lineWidth/2,$blue);
	$image->filledRectangle($label,$height-$threshold*$timesY-$lineWidth/2,$label+5*$lineWidth/2,$height-$threshold*$timesY+$lineWidth/2,$blue);
	$label+=4*$lineWidth;
}

#$threshold=6;
#$label=$extend;
#while($label<$width){
#	$image->rectangle($label,$height-$threshold*$timesY-$lineWidth/2,$label+5*$lineWidth/2,$height-$threshold*$timesY+$lineWidth/2,$red);
#	$image->filledRectangle($label,$height-$threshold*$timesY-$lineWidth/2,$label+5*$lineWidth/2,$height-$threshold*$timesY+$lineWidth/2,$red);
#	$label+=4*$lineWidth;
#}
#
#$threshold=5;
#$label=$extend;
#while($label<$width){
#	$image->rectangle($label,$height-$threshold*$timesY-$lineWidth/2,$label+5*$lineWidth/2,$height-$threshold*$timesY+$lineWidth/2,$yellow);
#	$image->filledRectangle($label,$height-$threshold*$timesY-$lineWidth/2,$label+5*$lineWidth/2,$height-$threshold*$timesY+$lineWidth/2,$yellow);
#	$label+=4*$lineWidth;
#}

$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",30,0,5,35,"$ARGV[2]");
binmode(IMAGE);
print IMAGE $image->png;
close IMAGE or die "Cant't close IMAGE: $!";
print "manhattan plot is finished...\n";









