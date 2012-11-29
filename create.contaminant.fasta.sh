#!/bin/bash
#===============================================================================
#
#          FILE:  create.contaminant.fasta.sh
# 
#   DESCRIPTION:  commands used to create the fasta contaminants file from /home/fgarcia/software/fastqc/FastQC/Contaminants/contaminant_list.txt
# 
#  REQUIREMENTS:  ---
#         NOTES:  ---
#        AUTHOR:  Fernando Garcia (), garcia@mpiib-berlin.mpg.de
#       COMPANY:  Molecular Biology at MPI
#       VERSION:  0.1
#       CREATED:  11/29/12
#      REVISION:  ---
#===============================================================================


# Removing the /r
cp /home/fgarcia/software/fastqc/FastQC/Contaminants/contaminant_list.txt /home/fgarcia/software/fastqc/FastQC/Contaminants/contaminant_list.txt.bck
dos2unix /home/fgarcia/software/fastqc/FastQC/Contaminants/contaminant_list.txt

# Creating fasta
cd /home/fgarcia/software/fastqc/FastQC/Contaminants/
grep -v "#" contaminant_list.txt | sed '/^$/d'| sed 's/ /-/g' | awk '{print ">" $1; print $2}' > /data/scratch/fastas/contaminants.fastqc.fa

# Creating bowtie index
cp /data/scratch/fastas/contaminants.fastqc.fa /data/scratch/ebwts/contaminants.fastqc.fa
cd /data/scratch/ebwts
bowtie-build contaminants.fastqc.fa contaminants.fastqc

# Creating bowtie2 index
cp /data/scratch/fastas/contaminants.fastqc.fa /data/scratch/ebwts/bowtie2/contaminants.fastqc.fa
cd /data/scratch/ebwts/bowtie2
bowtie2-build contaminants.fastqc.fa contaminants.fastqc