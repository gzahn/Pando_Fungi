#------------------------------------------------------------------#
# Download raw sequence data from SRA Project PRJNA1076591         #
#                                                                  #
# These are single-end PacBio reads                                #
# Requires that sra-toolkit be installed and in your system path   #
# (https://github.com/ncbi/sra-tools)                              #
#------------------------------------------------------------------#


# Setup ####
library(tidyverse)
fung <- read_csv("./Data/Combined_Metadata_2025.csv")


# Download fungal amplicons ####
accessions <- fung$sra_accession
filenames_f <- file.path(str_remove(fung$fwd_filepath,".gz"))
dl_filenames_1 <- paste0(accessions,"_1.fastq")


# run download in a for-loop for all accessions
# compress and rename and move files
if(!dir.exists("./Data/Raw")){
  dir.create("./Data/Raw",recursive = TRUE)
}
i=1
for(i in seq_along(accessions)[1]){
  system2("fasterq-dump",args = accessions[i])
  system2("gzip", args = c(dl_filenames_1[i]))
  file.rename(paste0(dl_filenames_1[i],".gz"),fung$fwd_filepath[i])
}



