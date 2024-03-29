#########################################
# Get ATAC fragments from ArchR Project #
#########################################

# Load libraries
library("ArchR")
library("qs")

# Set datasets
datasets <- c("timecourse", "Neuron", "NPC", "cocultured28")

# get fragments for each dataset
# iterate over them and execute the extraction
for (dataset in datasets){
  
  print(dataset)

  # Set ArchR Project directory
  archr_dir <- paste0("/g/scb/zaugg/marttine/dariaMultiome/ArchR/", dataset, "/ArchRproject.qs")
  
  # Load ArchR Project
  ArchRProj <- qread(archr_dir)
  
  # Get the fragments
  fragments <- getFragmentsFromProject(
    ArchRProj = ArchRProj,
  )
  
  # Have a look at them
  fragments[[1]]
  
  # Save them into a tsv file
  write.table(x = fragments[[1]], file = paste0("/g/scb/zaugg/deuner/SCENIC+/inputdata/", dataset, "_fragments.tsv"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\t")
  
}

# FORMAT THEM
# Then convert the strand column from "*" to "+" by; sed 's/\*/+/g' timecourse_fragments.tsv > timecourse_fragments_translated.csv
# Reorder the columns by: awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5 }' timecourse_fragments_translated.csv > timecourse_fragments_translated_ordered.tsv
# Change "#" for "_": sed 's/\#/_/g' timecourse_fragments_translated_ordered.tsv > timecourse_fragments_translated_ordered2.tsv
# Sort them: sort -t$'\t' -k1,1 -k2,2n -k3,3n timecourse_fragments_translated_ordered2.tsv > timecourse_fragments_translated_ordered3.tsv
# BED format: http://www.ensembl.org/info/website/upload/bed.html
# Remove "timecourse_": awk '{ gsub(/timecourse_/,"", $4); print } ' timecourse_fragments_translated_ordered3.tsv > timecourse_fragments_translated_ordered4.tsv
# Set tab as column delimiter: sed 's/ /\t/g' timecourse_fragments_translated_ordered4.tsv > timecourse_fragments_translated_ordered5.tsv
