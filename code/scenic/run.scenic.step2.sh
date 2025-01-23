#!/bin/bash

#### note: lustre file system doesn't support flock.
#### go to GPFS
cDir=`pwd`
outDir=$cDir/OUT.scenic/grn
shDir=$cDir/sh.scenic
mkdir -p $shDir

ls $cDir/OUT.scenic/sce/*.loom \
    | perl -ane 'chomp; /.+scenic.(.+?).mini.(.+?).loom$/; print "$1.$2\t$_\n"; ' \
    > $cDir/loom.list

while read aid lfile
do
	(
	cat<<-HERE
#!/bin/bash
#SBATCH -p q.all
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o S.$aid.%j.out
#SBATCH -e S.$aid.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
$cDir/../run.pyScenic.sh \\
    $lfile \\
    $outDir/$aid
HERE
)>$shDir/pyScenic.$aid.sh
done<$cDir/loom.list

