use strict;
use warnings;

###perl program/stat.allSNPSignificantPos.inOneDir.pl EMMAx.Result.2nd/Luohe kin Sesame_705lines.Luohe.2nd.peakSNP Sesame_705lines
if(!-e 'sigSNP'){
	system "mkdir sigSNP";
}
my @file=`ls -l $ARGV[0]/$ARGV[3]_*_chr01.$ARGV[1]| awk '{print \$9}'`;
open OUTPUT,">sigSNP/$ARGV[2].list" or die "Can't open the OUTPUT file...\n";
print OUTPUT "Trait\tLinkageGroup\tPosition\tpValue\n";

sub findSNPSig{
	my($posA,$posB);
	my(@value,$height,$pos);
	my $extend=1500000;		##DT??¡êo¡¤???????¨®¨°
	$pos=`awk '{if(\$2~/^[0-9]/ && \$3>max)(pos=\$2)(max=\$3);} END{print pos}' $ARGV[2].temp`;
	chomp($pos);
	
	$posA=$pos-$extend;
	$posB=$pos+$extend;
	system "awk '{if(\$2~/^[0-9]/ && (\$2<$posA || \$2>$posB)) print \$0}' $ARGV[2].temp >$ARGV[2].temp_filter";
	system "cp $ARGV[2].temp_filter $ARGV[2].temp";
	
	$height=`awk '{if(\$2~/^[0-9]/ && \$3>max)(max=\$3);} END{print max}' $ARGV[2].temp`;
	chomp($height);
	return($height,$pos);
}

my $threshold=5;		##DT??¡êo?D?¦Ì
my $range="3mb";	##DT??¡êo¡¤???????¨®¨°
my (@allPosList,$maxNum,$height,@value,$posA,$posB,$num,@position,$max_value,$pos_value,);
foreach(@file){
	if($_=~/$ARGV[0]\/$ARGV[3]_(\S+)_chr01.($ARGV[1])$/){
		chomp;
		print "$_\n";
		my $chr=1;
		@allPosList=();
		$maxNum=0;
		while($chr<=12){
			if($chr<10){
				system "cp $ARGV[0]/$ARGV[3]_$1_chr0$chr.$2 $ARGV[2].temp";
			}else{
				system "cp $ARGV[0]/$ARGV[3]_$1_chr$chr.$2 $ARGV[2].temp";
			}			
			$height=`awk '{if(\$2~/^[0-9]/ && \$3>max)(max=\$3);} END{print max}' $ARGV[2].temp`;
			chomp($height);
			$num=0;
			@position=();
			loopFIND: while($height>=$threshold){
				my $MaxHeight=$height;
				($max_value,$pos_value)=&findSNPSig;
				$height=$max_value;
				$position[$num]=$pos_value;
				print OUTPUT "$1\t$chr\t$pos_value\t$MaxHeight\n";
				$num++;
				if($max_value eq ""){
					last;
				}
			}
			if(@position>$maxNum){
				$maxNum=@position;
			}
			$num=0;
			foreach(@position){
				$allPosList[$num][$chr-1]=$_;
				$num++;
			}
			$chr++;
		}
	}
}
close OUTPUT;
system "rm $ARGV[2].temp $ARGV[2].temp_filter";
