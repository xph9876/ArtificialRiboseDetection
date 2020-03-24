# ArtificialRiboseDetection
Detect the artificial rNMPs captured by _ribose-seq_ and Ribose-Map.

## Introduction
Some fragments after restriction enzyme digestion don't contain the rNMP. However, they still remains after the whole _ribose-seq_ and Ribose-Map protocol due the imperfect digestion of T5 nuclease. That makes up a part of background noise. Some of them may have an additional dAMP at the end because of the dA-tailing step. 
This software is designed for background noise detection and subtraction if needed. It counts the number of the rNMPs incorporated in restriction enzyme cut site with or without dA tailing, calculates the ratio of them and subtracts them from the BED file.

## Usage
__subtraction.py__ is the main script to run the subtraction.

Inputs:
1. __list__                  list of calculation, fs species RE1,RE2,....
2. __raw_reads__             folder of raw reads

Parameters:
1. __-o O__                  output file
1. __--subtracted_output__ SUBTRACTED_OUTPUT                     Folder to output subtracted files, default=same to bed
1. __--ns__                  Do not do subtraction
1.  __--bed BED__             Folder of bed files, default=current folder
1.  __--mitochondrial_name MITOCHONDRIAL_NAME__
                        Name of mitochondrial DNA in genome file, default=chrM
1.  __--RElist RELIST__       pool of RElist, Name Pattern Pos, default=res_all.list
1.  __--genome GENOME__       Folder of genome, provide if needed
1.  __--threads THREADS__     Number of threads

