## replication timing data for 300 individuals at a variety of genomic windows specified in hg19 coordinates
https://www.thekorenlab.org/data > "DNA replication timing of 300 iPSC lines (HipSci)"

## lifting-over from hg19 to hg38 coordinates
I used the UCSC's CLI to perform the lifover, as the file is > 500MB in size.
The steps to run the CLI are described here: https://genome.ucsc.edu/FAQ/FAQdownloads.html#liftOver 
However, Brent has automated all the steps here: https://gist.github.com/brentp/894555/f23d1d6e0c988d6711acf2fe1a5bb930c3a19604
I used Brent's approach by defining an alias "lift", as described in his gist. 

## data location 
/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/replication-timing
