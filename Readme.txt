How to run CloverInJulia:

There are three files that need to be provided as input for the program.
The input sequence in fasta format
The motif library. An example is included at the end of this document.
The background sequence in fasta format.

Put those files inside the same folder as the julia files and edit the globalVairable.jl file. Change the file name global variables to the ones you provided. 

In addition, there are four other parameters that can be tuned to adjust the performance of CloverInJulia:
Hit_thresh: the lower bound for the raw_score of a result to be counted as a hit.
PseudoCount: used when normalizing the motif matrix to change it from a frequencies matrix to a probability matrix. 
Pthresh: the lower bound for the p value of a result for it to be counted as significant and outputted as the final output.
Shuffle: defines the number of shuffles that are conducted, this will affect the motif and p-value you get.

How to interpret the results:
A summary will be generated as a csv file recording all the significant motifs identified in the sequence along with their pvalues. For each sequence, a CSV file will be generated recording all motifs identified in the sequence along with the location of the motif and the nucleotides in the sequence at the location. Command line output are the same as the csv output except that in command line you will also be able to see a restatement of the input variable used. And if there are too many motifs identified, the results table might be abbreviated.

small changes
