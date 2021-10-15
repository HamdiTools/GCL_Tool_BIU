# GCL_Tool_BIU
Calculate the GCL of the input data with bootstrap

# Requirements:
1) python 3
2) Packages for python (if you don't have one of the following, right click inside python on the missing import -> Show Context Actions -> Install Package):<br /> numpy, math, random, scipy.

# How to include in your own python project:
1) put the gcl_library.py file insied your python project.
2) add the import: import gcl_library as gcl_lib
3) for example adding a bootsratp command to get a vector, simply write: gcl_lib.bootstrap(the required arguments)

# Usual running time:
For bootsrtap process: num_divisions=50 and bootsrtraps =100 on the two given files in this reposetory: around 8 minutes.
<br />For bootsrtap process: num_divisions=70 and bootsrtraps =70 on the two given files in this reposetory: around 9 minutes.
<br />For regular gcl value: num_divisions=100 on the two given files in this reposetory:around 20 seconds.

# How To run the example file :
download all the files from this github reposatory to a floder on your computer.<br />
A) Using Python Console, run the program with the following arguments:
1) data_arr [mandatory] - files array data with all the files in the current folder to make the GCL on.
2) num_division [optional - default:50] - Number of random gene division iterations for calculation.
3) bootsrtaps [optional - default:100] - Number of boot straps iterations for calculation.
4) bootstrap_percentage [optional - default:0.8] - Percentage of cells to choose for bootstrapping, a number between 0 and 1.
5) task_option [optional -default:'bootsrtap'] - either 'bootsrtap' or 'regular_calc' for the requested task.
* for default values - enter 'defuault' instead of the requested value.<br />An input for example: young.csv,old.csv 50 100 default bootstrap
