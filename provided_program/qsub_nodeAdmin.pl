use strict;
use warnings;

if(@ARGV!=2){
	print "Error:		perl program/qsub_nodeAdmin.pl emmax.sh \"command line\"  \n";
}
open SH,">$ARGV[0]" or die "Can't open the SH file: $ARGV[1] ...\n";
print SH "#!/bin/bash\n";
print SH "#PBS -N $ARGV[0]\n";
print SH "#PBS -l nodes=1:ppn=1\n";
print SH "#PBS -j oe\n";
print SH "#Run your executable\n";
print SH "echo starting\n";
print SH "date\n";
print SH "cd \$PBS_O_WORKDIR\n";
print SH "###################\n";
print SH "$ARGV[1]";
print SH "\n";
print SH "###################\n";
print SH "date\n";
print SH "echo ending\n";
close SH;

my $sub="qsub -l mem=1gb $ARGV[0]";
print "$sub\n";
system $sub;
sleep 1;
