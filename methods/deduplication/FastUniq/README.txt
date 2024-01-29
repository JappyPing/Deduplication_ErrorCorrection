Introduction:
=========================================================================
FastUniq as an ultrafast de novo tool for removal of duplicates in paired
short DNA sequence reads in FASTQ format. FastUniq identifies duplicates
by comparing sequences between read pairs and does not require complete
genome sequences as prerequisites. FastUniq is capable of simultaneously
handling reads with different lengths and results in highly efficient 
running time.


Installation:
=========================================================================

1). Make sure the gcc compiler installed on your computer (Version 4.0 or 
    above is recommanded).
   
2). Download the latest source code package of FastUniq
    (e.g. FastUniq-1.1.tar.gz),  and uncompress this package. 
   
3). Open terminal window, and go to "source" folder of the FastUniq. Open
    the "makefile" file, go to the "GCC_OPTION" line which is used to
    define the compiler arguments. Your can alter it following the gcc 
    compiler option's instructions as you needed. 

4). Type "make"

5). Now, "fastuniq" located in "source" folder is ready to use, you can
    move it to any location as you need. 


Unistall:
==========================================================================

To uninstall FastUniq, remove the "fastuniq" file located in the "source"
folder, or the "fastuniq" file moved to any location. 


FastUniq Program parameters:
==========================================================================

-i : The input file list of paired FSATQ sequence files [FILE IN]
        Maximum 1000 pairs

     This parameter is used to specify a list of paired sequence files in 
     FASTQ format as input, in which two adjacent files with reads in the
     same order belong to a pair.

-t : Output sequence format [q/f/p]
        q : FASTQ format into TWO output files
        f : FASTA format into TWO output files
        p : FASTA format into ONE output file
        default = q

     This parameter is used to specify sequence format in output file(s).
     FastUniq could output read pairs into two files in either FASTQ [q]
     or FASTA [f] format, in which reads in the same order belonging to a 
     pair. FastUniq could also output read pairs into a single file in 
     FASTA format [p], in which adjacent reads belonging to a pair.
     
-o : The first output file [FILE OUT]

-p : The second output file [FILE OUT]
     Optional. ONLY required when output sequence format(-t) is specify as
     [q] or [f]. 

-c : Types of sequence descriptions for output [0/1]
        0 : The raw descriptions
        1 : New serial numbers assigned by FastUniq
        default = 0
