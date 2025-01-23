#!/bin/bash

Rscript ./w.CD8.R
Rscript ./w.CD4.R

### submit the jobs

###

ls OUT.gsea/C2.CP/gsea.*.classic.GSEA/*.GseaPreranked*/gsea_report_for_na_*.xls \
	| perl -ane 'chomp;/\/gsea.(.+?)\.(.+?).classic.GSEA\//; print "$1\t$2\t$_\n" ' \
	> gsea.tb.C2.CP.list

ls OUT.gsea/KEGG.flt/gsea.*.classic.GSEA/*.GseaPreranked*/gsea_report_for_na_*.xls \
	| perl -ane 'chomp;/\/gsea.(.+?)\.(.+?).classic.GSEA\//; print "$1\t$2\t$_\n" ' \
	> gsea.tb.KEGG.flt.list

(
while read stype mcls afile
do
	awk -F"\t" -v OFS="\t" 'NR==1{print "stype","meta.cluster",$0} NR>1{print "'$stype'","'$mcls'",$0}' $afile 
done<gsea.tb.C2.CP.list
) | awk 'NR==1||!/^stype/' >gsea.tb.C2.CP.txt

(
while read stype mcls afile
do
	awk -F"\t" -v OFS="\t" 'NR==1{print "stype","meta.cluster",$0} NR>1{print "'$stype'","'$mcls'",$0}' $afile 
done<gsea.tb.KEGG.flt.list
) | awk 'NR==1||!/^stype/' >gsea.tb.KEGG.flt.txt


