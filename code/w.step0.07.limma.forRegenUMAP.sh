#!/bin/bash


cDir=`pwd`
sDir=`R --slave -e 'sDir <- system.file("script",package="scPip"); cat(sDir)'`
echo $sDir

for group_mode in "multiAsTwo"
do
        for ss in CD8
        do
        	while read dataID measurement platform seufile scefile
        	do
                aid=$dataID
				a_res="meta.cluster"
				#a_res="meta.cluster.coarse"
				shDir=./sh.limma.inte.post.$ss.$a_res.$group_mode.continumOnly.v1
				mkdir -p $shDir

        	(
        	cat <<-HERE
#!/bin/bash
#SBATCH -p q.all
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-node=4
#SBATCH -o S.$dataID.%j.out
#SBATCH -e S.$dataID.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
$sDir/wrapper.run.limma.R \\
	-b $scefile \\
	-o $cDir/OUT_Fig2/regenUMAP/limma.post.$a_res.$group_mode.continumOnly.v1/T.$ss.$aid/T.$ss.$aid \\
	--platform $platform \\
	--group $a_res \\
	--filter CD8.c03.Tm.RPS12,CD8.c04.Tm.CD52,CD8.c17.Tm.NME1,CD8.c07.Temra.CX3CR1,CD8.c08.Tk.TYROBP,CD8.c09.Tk.KIR2DL4,CD8.c15.ISG.IFIT1,CD8.c16.MAIT.SLC4A10,CD8.c13.Tex.myl12a,CD8.c14.Tex.TCF7 \\
	--groupMode $group_mode \\
	-n 6 \\
	-m $measurement
HERE
	    	)>$shDir/limma.$ss.$aid.$a_res.$group_mode.sh
		done< <(awk '!/^#/ && !/^data.id/' list/sce.CD8.fullPath.list)
	    done
done

##### submit the jobs
# cd $shDir
# ls *.sh | awk '{print "sbatch "$0}' | bash

