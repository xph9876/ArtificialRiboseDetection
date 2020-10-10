# ArtificialRiboseDetection
Detect the artificial rNMPs captured by _ribose-seq_ and Ribose-Map.

## Introduction
Some fragments after restriction enzyme digestion don't contain the rNMP. However, they still remains after the whole _ribose-seq_ and Ribose-Map protocol due the imperfect digestion of T5 nuclease. That makes up a part of background noise. Some of them may have an additional dAMP at the end because of the dA-tailing step. 
This software is designed for background noise detection and subtraction if needed. It counts the number of the rNMPs incorporated in restriction enzyme cut site with or without dA tailing, calculates the ratio of them and subtracts them from the BED file.

## Citation
Balachander, S., Gombolay, A. L., Yang, T., Xu, P., Newnam, G., Keskin, H., El-Sayed, W., Bryksin, A. V., Tao, S., Bowen, N. E., Schinazi, R. F., Kim, B., Koh, K. D., Vannberg, F. O., & Storici, F. (2020). Ribonucleotide incorporation in yeast genomic DNA shows preference for cytosine and guanosine preceded by deoxyadenosine. _Nature communications, 11_(1), 2447. https://doi.org/10.1038/s41467-020-16152-5

## Dependency
Following softwares are needed:
- bedtools
- Linux build-in __grep__ and __wc__

This software only requires Python3 standard libraries.

## Usage
__subtraction.py__ is the main script to run the subtraction.

Inputs:
1. __list__                  List of files whose background noise ratio should be calculated. This list contains 3 columns:
   1. Name: Name of the library. Should be same with BED file and FASTQ file.
   2. Reference genome: Name of the reference genome. Corresponding FASTA file and FASTA index file (build by samtools) are needed in the __genome__ folder.
   3. Restriction enzyme used. If multiple restrction enzymes are used, they should be separated by "," with no space.
   Please find __test.list__ as an example.
2. __raw_reads__             Folder of FASTQ files. The extension name should be __".fq"__ rather than ".fastq".

Parameters:
1. __-o O__                  Output file name.
1. __--subtracted_output SUBTRACTED_OUTPUT__                     Folder to output subtracted BED files, default = same to __bed__
1. __--ns__                  Do not do subtraction, calculation only.
1.  __--bed BED__             Folder of BED files, default = current folder.
1.  __--special S [ S ...] __
                        Name of special chromosomes in genome file, which is calculated separately and excluded from nuclear DNA. default = [chrM]
1.  __--RElist RELIST__       Pool of RElist, default = res_all.list. This list have 3 columns.
    1. Name of the restriction enzyme.
    1. Recognition pattern of the restriction enzyme. Degraded bases are supported.
    1. Cut after which base.
    Please find __res_all.list__ as an example.
1.  __--genome GENOME__       Folder of genome, provide if needed.
1.  __--threads THREADS__     Number of threads

## License
This software is under GNU GPL v3.0 license.

## Contact
If you have any question, please contact me at <pxu64@gatech.edu>.

