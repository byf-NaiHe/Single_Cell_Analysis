#!/bin/bash

cDir=`pwd`

binPath=$cDir/"wrapper.reClassify.R"

###for m_mode in "rf" "svm"

for m_mode in "svm"
do
        for ss in CD8 CD4
        do
        while read data.id measurement     platform        seufile scefile
	    do
	        a_res="meta.cluster"
		shDir=./sh.reClassify.v1.$ss.$m_mode
		mkdir -p $shDir

        	(
        	cat <<-HERE
#!/bin/bash
#SBATCH -p q.all
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o S.$ss.$m_mode.%j.out
#SBATCH -e S.$ss.$m_mode.%j.err
#SBATCH --no-requeue
echo \`hostname\`
$binPath \\
	-i $scefile \\
    -o $cDir/../data/expression/$ss/integration/reClassify.v1/T.$ss.$aid/T.$ss.prob.mat.svm.$aid.rds \\
	--id $aid \\
	--mode $m_mode \\
	--ncore 8
HERE
	    	)>$shDir/reClassify.$ss.$aid.$m_mode.sh
        done< <(awk 'NR>1 && !/^#/' $cDir/list/sce.$ss.fullPath.list)
	    done
done

########## submit jobs to run the scripts ########
#cd $shDir
#ls *.sh | awk '{print "sbatch "$0}' | bash
#cd $cDir
##################################################




