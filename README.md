# Deduper

Given a sorted SAM file of uniquely mapped reads, this tool will remove all PCR duplicates and retain only a single copy of each read. Before the tool is used the SAM file must be sorted. This tool accounts for:
- all possible CIGAR strings and adjusts for soft clipping
- Strand
- Single-end reads
- Known UMIs

In the future, I will add functionality to handle:
    - Single-end vs paired-end reads
    - Known UMIs vs randomers
    
The output of the tool is an a SAM file with "_deduped" appended to the filemane. The file will contains the first read encountered if duplicates are found

