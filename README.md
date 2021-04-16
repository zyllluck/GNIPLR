#GNIPLR

**GNIPLR** (gene networks inference based on projection and lagged regression) is a method for constructing  gene regulatory network.
Different data needs to train different parameters, the gniplr.py only provides one data parameter. If you need to run the methods, you can contact the author to get other parameters.


## Run GNIPLR

gniplr.py


## Input Data Format
The example 'BLCA_cancer.txt' contains 10 genes, 414 samples.
The text 'BLCA_cancer_gold' is the network structure.

## Output file format
Four-column matrix file, each column is the regulator, target gene, threshold, and gold standard edge


