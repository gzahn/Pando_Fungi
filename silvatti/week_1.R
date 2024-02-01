library(dada2)
library(ShortRead)
library(Biostrings)
library(tidyverse)

#Cmd shift m %>% 
# option - <-

fq <- ShortRead::readFastq("./data/2554_pass_1.fastq.gz")
fq@sread
fq@quality
fq@id
ShortRead::reverseComplement(fq) %>% 
  ShortRead()