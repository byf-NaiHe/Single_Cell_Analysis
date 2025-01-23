#!/bin/bash

#Rscript scenic.adj.comb.R

cDir=`pwd`
outDir="$cDir/OUT.scenic/grn.comb"
shDir=`pwd`/sh.scenic.comb
mkdir -p $shDir

############# step 3 ##############
Rscript run.scenic.step3.comb.adj.R

cd $cDir/OUT.scenic/grn/
ln -s comb.CD4.adj.csv.gz comb.zhangLab10X.CD4.adj.csv.gz
ln -s comb.CD4.adj.csv.gz comb.zhangLabSS2.CD4.adj.csv.gz
ln -s comb.CD8.adj.csv.gz comb.zhangLab10X.CD8.adj.csv.gz
ln -s comb.CD8.adj.csv.gz comb.zhangLabSS2.CD8.adj.csv.gz
cd $cDir

while read aid lfile
do
	(
	cat<<-HERE
#!/bin/bash
#SBATCH -p fat4way
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -o S.$aid.%j.out
#SBATCH -e S.$aid.%j.err
#SBATCH --no-requeue
#SBATCH -A zeminz_g1
#SBATCH --qos=zeminzf4w
source /lustre1/zeminz_pkuhpc/zhenglt/.bashrc_zhenglt
echo \`hostname\`
$cDir/../run.pyScenic.sh \\
	-a $cDir/OUT.scenic/grn/comb.$aid.adj.csv.gz \\
    $lfile \\
    $outDir/$aid
HERE
)>$shDir/pyScenic.$aid.sh
done<$cDir/loom.list

####
## submit the job
###

############# step 4 ##############
### some python
python run.scenic.step4.zhangLab10X.CD8.vComb.py
python run.scenic.step4.zhangLab10X.CD4.vComb.py
python run.scenic.step4.zhangLabSS2.CD8.vComb.py
python run.scenic.step4.zhangLabSS2.CD4.vComb.py

### some R
Rscript run.scenic.step6.comb.R
#Rscript run.scenic.step5.zhangLabSS2.vComb.R
#Rscript run.scenic.step5.zhangLab10X.vComb.R



