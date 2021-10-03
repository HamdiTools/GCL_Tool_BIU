# GCL_Tool_BIU
Calculate the GCL of the input data with bootstrap

# How To:
A) Using Python Console, run the program with the following arguments:
[1] mandatory: a list of csv files with commas between each file.
[2] optional: number of random gene division/ bootsrtap iterations for calculation. (default will be 100)
[3] optional: percentage of cells to choose for bootstrapping. (default will be 0.8).
[4] optional: either 'bootsrtap' or 'regular_GCL' for the requested task. (default will be 'bootsrtap').
* for default values - enter defuault instead of the requested value, an input for example:
young.csv,old.csv 100 default bootstrap

B) using jupiter notebook:
There are two examples in the file - the first for two files bootsrtap and the latter foor a single regular_gcl calculation
