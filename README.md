# CSER
geneexp.txt: Expression levels of genes in samples, with row and column names removed. Used as an input file for MATLAB.

CSER1: Code for calculating causal strength.

MATLAB output:

net1.xlsx: Indicates whether there is an interaction between genes, where 0 represents no interaction and 1 represents interaction.
netvalue1.xlsx: Values of causal strength between genes.
If a gene is not related to any other gene, it should be removed.

CSER2.py: Code for constructing gene regulatory networks.

intersectgeneEXP.txt: Input file for Python, containing expression levels of genes in samples after removing independent genes.

CC_result_network.txt: Output file for the gene regulatory network, where the first column is the gene, the second column is the target gene, and the third column is the weight (which can be positive or negative).
