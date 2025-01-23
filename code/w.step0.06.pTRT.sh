#!/bin/bash

cDir=`pwd`

################ CD8+ ###################################
shDir=./sh.gsea.perCell.zscore.mini.CD8
mkdir -p $shDir
while read dataID measurement platform seufile scefile
do
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
$cDir/wrapper.pTRTs.step1.gsea.perCell.R \\
    -i $scefile \\
    -o $cDir/OUT_Fig1/pTRT/pTRT.signaling.CD8 \\
    -c $cDir \\
    -d $dataID
HERE
)>$shDir/gsea.$dataID.sh
done< <(sed '1,1d' list/sce.CD8.fullPath.list)

################ CD4+ ###################################
shDir=./sh.gsea.perCell.zscore.mini.CD4
mkdir -p $shDir
while read dataID measurement platform seufile scefile
do
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
$cDir/wrapper.pTRTs.step1.gsea.perCell.R \\
    -i $scefile \\
    -o $cDir/OUT_Fig1/pTRT/pTRT.signaling.CD4 \\
    -c $cDir \\
    -d $dataID
HERE
)>$shDir/gsea.$dataID.sh
done< <(sed '1,1d' list/sce.CD4.fullPath.list)


