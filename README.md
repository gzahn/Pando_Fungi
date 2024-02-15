# Foliar Fungi of the Pando Aspen Clone (USA: Utah)


**Foliar endophytic and epiphytic fungi were sampled from healthy trees across the Pando aspen clone**

### Raw data available on the Sequence Read Archive under accession [PRJNA1076591](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1076591)






___


Samples were prepared for ITS2 amplicon sequencing using a modified version of the Illumina 16S Metagenomics Sequencing Protocol: https://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf

Briefly:
gDNA samples were amplified using Q5 High-Fidelity 2X Master Mix (NEB) and the following primers, consisting of the ITS3F and ITS4R binding sites plus overhangs compatible with Illumina Unique Dual Indexes:

5' TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGCATCGATGAAGAACGCAGC 3' (GCATCGATGAAGAACGCAGC = ITS3F binding sequence)

5' GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTCCTCCGCTTATTGATATGC 3' (TCCTCCGCTTATTGATATGC = ITS4R binding sequence)

Amplification PCR conditions:
98°C	30sec	
98°C	10sec	
61°C	30sec	
72°C	25sec	
(25 cycles of three previous steps)
72°C	2min	
10°C	hold

PCR reactions were purified using SPRI bead cleanup, then used as input for the indexing PCR reaction. Samples were indexed using Q5 High-Fidelity 2X Master Mix (NEB) and IDT for Illumina DNA/RNA UD Indexes, Tagmentation (Illumina).

Indexing PCR conditions:
98°C	30sec
98°C	10sec
56°C	30sec
72°C	25sec
(8 cycles of the three previous steps)
72°C	2min
10°C	hold

PCR reactions were purified using SPRI bead cleanup, then pooled and diluted for sequencing.

Sequencing was performed on the Illumina NextSeq2000 platform using a 600 cycle flow cell kit to produce 2x300bp paired reads. 30-40% PhiX control (unindexed) was spiked into the library pool to support optimal base calling of low diversity libraries on patterned flow cells.

Read demultiplexing, read trimming, and run analytics were performed using DRAGEN v3.10.12, an on-board analysis software on the NextSeq2000. We include fastqc metrics as a best practice and for examination in the case of unexpected outputs.

