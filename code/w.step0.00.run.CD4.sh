#!/bin/bash

g_prj_id="panC"

cDir=`pwd`
### columns required in $file_obj_in: data.id measurement platform seufile scefile
### columns required in seurat or sce object: "patient"
file_obj_in="$cDir/list/sce.CD4.list"
file_obj_in_fullPath="$cDir/list/sce.CD4.fullPath.list"
file_int_in="$cDir/list/sce.CD4.int.list"
file_limma_in="$cDir/list/sce.CD4.forLimma.list"
file_limma_out="$cDir/list/sce.CD4.limma.sc.list"
cellType="CD4"
oDir="$cDir/OUT.byDataset"
shDir="./sh.byDataset.$cellType"
mkdir -p $shDir

sed 's#../../#'$cDir/../'#' $file_obj_in > $file_obj_in_fullPath

sDir=`R --slave -e 'sDir <- system.file("script",package="scPip"); cat(sDir)'`
echo $sDir

gene_int_file="$cDir/../data/expression/CD4/integration/int.CD4.S35.gene.G1500.tb"

######### 1.1 generate scripts to run Seurat for each dataset
while read aid measurement platform seufile scefile resolution
do
(
cat <<-HERE
#!/bin/bash
#SBATCH -p q.all
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o S.$cellType.$aid.%j.out
#SBATCH -e S.$cellType.$aid.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
$sDir/run.seurat3.basic.R \\
	-a $seufile \\
	-b $scefile \\
	-o $oDir/$cellType/$cellType.$aid/$cellType.$aid \\
	-d 15 \\
	-n 12 \\
	--resolution 2 \\
    --deg \\
	-m $measurement \\
	--platform $platform
HERE
)>$shDir/byDataset.$cellType.$aid.sh
done< <(awk 'NR>1' $file_obj_in_fullPath)

########## 1.2 submit jobs to run the scripts ########
#cd $shDir
#ls *.sh | awk '{print "sbatch "$0}' | bash
#cd $cDir
##################################################

######## 2.1 prepare file list for integration ############
sed '1,1d' $file_obj_in_fullPath \
    | perl -ane 'BEGIN{ print "data.id\tmeasurement\tplatform\tdefile\tscefile\tseufile\n"  }
                chomp; print join("\t",@F[0..2],
                "'$oDir'/'$cellType'/'$cellType'.$F[0]/limma/'$cellType'.$F[0].de.out.limma.rda",
                "'$oDir'/'$cellType'/'$cellType'.$F[0]/'$cellType'.$F[0].sce.rds",
                "'$oDir'/'$cellType'/'$cellType'.$F[0]/'$cellType'.$F[0].seu.rds")."\n" ' \
    > $file_int_in
########################################################


######################## 2.2 first run ########################
#### $gene_int_file contains genes used to obtain ../data/expression/CD4/integration/int.CD4.S35.sce.merged.rds
#### wrapper.run.inte.R with 
####   --geneFile $gene_int_file \\
#### wii yield result very similar with original result, but not exact the same
#### because of differences of the computer enviroment, such as R version, Seurat version etc.

(
cat <<-HERE
#!/bin/bash
#SBATCH -p q.all
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o S.int.$cellType.%j.out
#SBATCH -e S.int.$cellType.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
$sDir/wrapper.run.inte.R \\
    --inFile $file_int_in \\
    --outPrefix $cDir/OUT.int.$cellType/int.$cellType \\
    --corVar S.Score,G2M.Score,DIG.Score1
HERE
)>$shDir/inte.$cellType.sh
########################

########## 2.3 submit jobs to run the scripts ########
#cd $shDir
#sbatch inte.$cellType.sh
#cd $cDir
##################################################

######################## 3. cluster annotation ########################
./w.step0.1.ann.$cellType.R

######################## 4.1 run limma per dataset ########################
join -1 1 -2 1 \
    <(cut -f 1-3 $file_int_in|sort -k 1r,1) \
    <(ls $cDir/OUT.int.$cellType/sce/*.sce.rds | perl -ane 'BEGIN{print "data.id\tscefile\n" } chomp; /sce\/(.+?).sce.rds/;print "$1\t$_\n"' | sort -k 1r,1) \
    | sed 's/\s\+/\t/g' \
    | awk -F"\t" -v OFS="\t" 'BEGIN{ print("data.id\tmeasurement\tplatform\tscefile") } !/^data.id/{print $0}' \
    > $file_limma_in

while read data_id measurement platform scefile
do
(
cat <<-HERE
#!/bin/bash
#SBATCH -p q.all
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -o S.limma.$cellType.$data_id.%j.out
#SBATCH -e S.limma.$cellType.$data_id.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
echo \`hostname\`
$sDir/wrapper.run.limma.R \\
        -b $scefile \\
        -o $cDir/OUT.int.$cellType/limma.sc/$data_id/limma.sc.$data_id \\
        --platform $platform \\
        --group "meta.cluster" \\
        --groupMode "multiAsTwo" \\
        -n 8 \\
        -m $measurement
HERE
)>$shDir/limma.$cellType.$data_id.sh
done < <(awk '!/^data.id/' $file_limma_in)

########## 4.2 submit jobs to run the scripts ########
#cd $shDir
# ls limma.$cellType.*.sh | awk '{print "sbatch "$0}' | bash
#cd $cDir
##################################################

######################## 5. limma to sce ########################
### please check that the format of data.id is: ^(cancerType).(dataset)$, and ther are no "." (dots) in cancerType and dataset 
join -1 1 -2 1 \
    <(cut -f 1-3 $file_int_in|sort -k 1r,1) \
    <(ls $cDir/OUT.int.$cellType/limma.sc/*/*.de.out.rda | perl -ane 'BEGIN{print "data.id\tdfile\n" } chomp; /.+limma.sc.(.+?).de.out/;print "$1\t$_\n"' | sort -k 1r,1) \
    | sed 's/\s\+/\t/g' \
    > $file_limma_out

$sDir/wrapper.convertLimmaToSCE.R \
    --limmaFile $file_limma_out \
    --outPrefix $cDir/OUT.int.$cellType/int.$cellType.limma.sc \
    --ncores 8

######################## 6. prepare data for web ########################
mkdir $cDir/OUT.data.web
(
cat <<-HERE
$cDir/OUT.int.$cellType/int.$cellType.meta.tb.rds 
$cDir/OUT.int.$cellType/int.$cellType.sce.merged.rds
$cDir/OUT.int.$cellType/int.$cellType.limma.sc.sce.pb.rds
$cDir/OUT.int.$cellType/int.$cellType.limma.sc.gene.desc.tb.rds
$cDir/OUT.int.$cellType/int.$cellType.limma.sc.geneTableLong.rds
$cDir/OUT.int.$cellType/int.$cellType.limma.sc.geneTableLong.collapsed.rds 
$cDir/OUT.int.$cellType/int.$cellType.colSet.rds
HERE
) | perl -ane 'chomp; ($dname,$bname)=/^(.+)\/(.+)$/; print "ln -s $_ '$cDir'/OUT.data.web/'$g_prj_id'.$bname\n"' \
  | bash

#,DatasetName,DatasetSource,perMiniCluster,perMetaCluster,meta.perCell,geneTableLong,geneDesc,colSet
printf "$g_prj_id.$cellType,$g_prj_id.$cellType,$cellType,$g_prj_id.int.$cellType.sce.merged.rds,$g_prj_id.int.$cellType.limma.sc.sce.pb.rds,$g_prj_id.int.$cellType.meta.tb.rds,$g_prj_id.int.$cellType.limma.sc.geneTableLong.collapsed.rds,$g_prj_id.int.$cellType.limma.sc.gene.desc.tb.rds,$g_prj_id.int.$cellType.colSet.rds\n" > $cDir/OUT.data.web/dataset_map.csv

############### final results ###################
# $cDir/OUT.int.$cellType/int.$cellType.limma.sc.gene.desc.tb.rds
# $cDir/OUT.int.$cellType/int.$cellType.limma.sc.sce.pb.rds
# $cDir/OUT.int.$cellType/int.$cellType.meta.tb.rds
# $cDir/OUT.int.$cellType/int.$cellType.seu.merged.rds
# $cDir/OUT.int.$cellType/int.$cellType.sce.merged.rds

### important columns: geneID median.F.rank
# OUT.int.$cellType/int.$cellType.gene.rank.tb.flt.rds



