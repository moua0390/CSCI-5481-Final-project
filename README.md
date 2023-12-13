# CSCI 5481 Final project (Group 22)

## Group members
+ Hmong Moua (`moua0390`)
+ Parsa Moradi (`moradi`)

## Proposal
We will attempt to perform Covid genome assembly (project option 2). This will involve applying quality trimming, correcting sequencing errors, dividing into k-mers (testing various values for k), and constructing as well as plotting De Bruijn graph on raw shotgun sequencing data of the SARS-CoV2 spike protein. Afterwards, we will attempt to assemble the genome by finding an Euclearian path, handling arising issues then align an assembled contig to the spike protein from the separate-genes file in Homework 1 and create a scatter plot to compare the location of their bases. Finally we will compare our assembly performance via known metrics (such as N50) with benchmark methods.
