#!/bin/bash

# set batch mode outputs directory
path="/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/"
cd $path

# takes as arguments dataset ($1) and correlation method ($2)

# create folder to store all the QC plots
mkdir $path${1}_batch_mode_${2}/peak_Gene_diagnosticPlots_all

# define cluster resolutions tested
res=("0.1" "0.5" "0.25" "0.75" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "12" "14" "16" "18" "20")

# se titerable variable
i=0

# iterate over folders
for dir in $path${1}"_batch_mode_"${2}"/Batch_Mode_Outputs/"
do
	cd $dir
	# iterate over cluster resolutions folders
	for dir2 in $(ls -v)
	do 
		cd ${dir2}"/plots/"
		cp peakGene_diagnosticPlots_all.pdf $path${1}_batch_mode_${2}/peak_Gene_diagnosticPlots_all/peak_Gene_diagnosticPlots_all${res[i]}.pdf
	#	pdftk /g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/$1_batch_mode_$2/peak_Gene_diagnosticPlots_all_MERGED.pdf peakGene_diagnosticPlots_all.pdf cat 1 output /g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/$1_batch_mode_$2/peak_Gene_diagnosticPlots_all_MERGED.pdf
		echo $dir2 "copied"
		cd $path${1}"_batch_mode_"${2}"/Batch_Mode_Outputs/"
		i=$i+1
	done
done

# merge all PDFs
cd $path${1}_batch_mode_${2}/peak_Gene_diagnosticPlots_all/
gs -dFirstpage=1 -dLastPage=1 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/$1_batch_mode_$2/peak_Gene_diagnosticPlots_all/peak_Gene_diagnosticPlots_all_MERGED.pdf $(ls -v)
