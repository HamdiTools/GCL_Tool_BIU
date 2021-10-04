# GCL_Tool_BIU
Calculate the GCL of the input data with bootstrap

# SetUp:
1) download all the files from this github reposatory to a floder on your computer.
2) make sure you have python 3.9 installed.

# How To:
A) Using Python Console, run the program with the following arguments:
<br />[1] mandatory: a list of csv files with commas between each file.
<br />[2] optional: number of random gene division/ bootsrtap iterations for calculation. (default will be 100)
<br />[3] optional: percentage of cells to choose for bootstrapping. (default will be 0.8).
<br />[4] optional: either 'bootsrtap' or 'regular_GCL' for the requested task. (default will be 'bootsrtap').
* for default values - enter defuault instead of the requested value, an input for example:
young.csv,old.csv 100 default bootstrap
<br />
B) using jupiter notebook:
<br />There are two examples in the file - the first for two files bootsrtap and the latter foor a single regular_gcl calculation.
<br />It is easier here to play with the arguments and get resaults according to the input you insert.
