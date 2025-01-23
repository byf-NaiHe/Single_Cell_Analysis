#!/bin/bash

optA=""

while getopts a: opt
do
	case $opt in 
	a)
		optA=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 <-a adj file> <loomfile> <out_prefix>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <-a adj file> <loomfile> <out_prefix>"
	exit 1
fi

loomfile=$1
out_prefix=$2

outDir=`dirname $out_prefix`
mkdir -p $outDir

echo begin at: `date`

cd $outDir

export HDF5_USE_FILE_LOCKING='FALSE'

ncores=8

SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")
echo $SCRIPT_DIR

if [ "$optA" != "" ] && [ -f $optA ];then
    echo "use previsouly saved adj.csv file"
	gzip -cd $optA > $out_prefix.adj.csv
else
    echo "to run: pyscenic grn"
	pyscenic grn \
		-o $out_prefix.adj.csv \
		-m grnboost2 \
		--seed 123456 \
		--num_workers $ncores \
		$loomfile \
        $SCRIPT_DIR/../data/external/scenic/transcription.factor.activity.GO0003700.symbol.list

		##--mask_dropouts \

fi

	##-a \

pyscenic ctx \
	-o $out_prefix.reg.csv \
	--expression_mtx_fname $loomfile \
    --num_workers $ncores \
	--annotations_fname $SCRIPT_DIR/../data/external/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
	$out_prefix.adj.csv \
	$SCRIPT_DIR/../data/external/scenic/hg19-500bp-upstream-7species.mc9nr.feather \
	$SCRIPT_DIR/../data/external/scenic/hg19-tss-centered-10kb-7species.mc9nr.feather

pyscenic aucell \
    $loomfile \
    $out_prefix.reg.csv \
    --output $out_prefix.pyScenic.loom \
	--seed 123456 \
    --num_workers $ncores

echo end at: `date`
