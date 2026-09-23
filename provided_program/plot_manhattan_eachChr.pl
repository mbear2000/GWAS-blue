#!/usr/bin/perl
##perl plot_manhattan_eachChr.pl <247line.4x.Genotype.chr> <_byR.GelConsistency> <output-filename> <text>
use warnings;
use strict;
use GD;
if(@ARGV!=3){
	print "\nError:  perl plot_manhattan_eachChr.pl 2nd3rd_20GrainWeight2012_Ave_chr01.hIBS 2nd3rd_20GrainWeight2012_Ave_hIBS 2nd3rd_20GrainWeight2012_Ave\n";
	print "!!!check reference file: /data9/home/yzhao/reference/IRGSP-1.0_genome.fasta.fai		\n\n";
	exit;
}

unless(-e "plot"){
	my $dir="mkdir plot";
	system $dir;
}

my $reference="/data9/home/yzhao/reference/IRGSP-1.0_genome.fasta.fai";
open REF,"<$reference" or die "Can't open the reference file: $reference\n\n";
my (@chrLen,@maxchrLen,@maxX,@value,$chr);
$maxX[0]=0;
while(<REF>){
	chomp;
	@value=split /\t/,$_;
	$chr=$value[0];
	$chr=~s/chr//;
	$chr=$chr+1-1;
	$chrLen[$chr]=$value[1];
	$maxchrLen[$chr]=$value[1]/1000000;
}
close REF;

#######################################
##ÐÞ¸Ä²ÎÊý
my $Xaxis=2000;
my $Yaxis=600;
my $radius=8;
my $extend=35;
my $fontSize=12;
my $lineWidth=3;
#######################################

my $threshold=5;
my $n=1;
while($n<=12){
	my $file_name=$ARGV[0];
	if($n<10){
		$file_name=~s/chr01/chr0$n/;
	}else{
		$file_name=~s/chr01/chr$n/;
	}
  open INPUT,"$file_name" or die "Can't open INPUT:$n $!";
	my @position=();
	my @pheno=();
  while(<INPUT>){
		chomp;
		if(!($_=~/^CHROM/)){
			@value=split /\s+/,$_;
			my $chr=$value[1]/1000000;
			push(@position,$chr);
			push(@pheno,$value[2]);
		}
	}
	close INPUT or die "Can't close INPUT: $!";
	
  my $maxX=$maxchrLen[$n];
  my @sort_pheno=sort{$b<=>$a}@pheno;
  my $maxY=1+int($sort_pheno[0]);  	
	my $timesX=$Xaxis/$maxX;
	my $timesY=$Yaxis/$maxY;
	
	my $outFile="plot/$ARGV[1]_chr01.png";
	if($n<10){
		$outFile=~s/chr01/chr0$n/;
	}else{
		$outFile=~s/chr01/chr$n/;
	}
  open IMAGE,">$outFile" or die "Can't open the image file :$n $!";
  my $image=GD::Image->new($Xaxis+$extend+$extend/2,$Yaxis+$extend*3);
  my $white = $image->colorAllocate(255,255,255);
  my $black = $image->colorAllocate(0,0,0 );
  my $red = $image->colorAllocate(255,0,0);
	my $width=$Xaxis+$extend;
	my $height=$Yaxis+$extend;
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

	my $num=0;
	my $length=@position;
	while($num<$length){
		my $x=$extend+$position[$num]*$timesX;
		my $y=$height-$pheno[$num]*$timesY;
		$image->arc($x,$y,$radius,$radius,0,360,$color{$n});
		$image->filledArc($x,$y,$radius,$radius,0,360,$color{$n});
		$num++;
	}
  
	$image->setThickness($lineWidth);
  ###	plot-threshold
	if($threshold<$maxY){
		$image->line($extend,$height-$threshold*$timesY,$width+$lineWidth,$height-$threshold*$timesY,$red);
	}	
  
	$image->line($extend,$height,$width+$lineWidth,$height,$black);
  $image->line($extend,$extend-$lineWidth,$extend,$height,$black);
  my $label=0;
  while($label<=$maxX){
  	if($label%5==0){
  		$image->line($extend+$label*$timesX,$height,$extend+$label*$timesX,$height+6*$lineWidth,$black);
	    $image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend+$label*$timesX-$fontSize/2,$height+8*$lineWidth+$fontSize,"$label");
  	}else{
    $image->line($extend+$label*$timesX,$height,$extend+$label*$timesX,$height+3*$lineWidth,$black);
    #$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend+$label*$timesX-$fontSize/2,$height+8*$lineWidth+$fontSize,"$label");
	}
    $label++;
  }
  $label=0;
  while($label<$maxY){
    $image->line($extend,$height-$label*$timesY,$extend-3*$lineWidth,$height-$label*$timesY,$black);
    $image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize,0,$extend-4*$lineWidth-$fontSize,$height-$label*$timesY+$fontSize/2,"$label");
    $label++;
  }
  if($n<10){
  	$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize*2,0,$width/2,$height+$extend+2*$fontSize+$lineWidth,"Chromosome 0$n");
  }else{
  	$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",$fontSize*2,0,$width/2,$height+$extend+2*$fontSize+$lineWidth,"Chromosome $n");
  }

  #$image->stringFT(color, font_path, size, angle, x, y, string)

  #$image->stringFT($black,"/data5/home/yzhao/program/times.ttf",2*$fontSize,0,$extend-$fontSize-$fontSize,$extend-$fontSize,"$ARGV[2].chr$n.png");
  binmode(IMAGE);
  print IMAGE $image->png;
  close IMAGE or die "Can't close IMAGE:$n $!";
  $n++;
}

print "each chromosome is finished...\n";
