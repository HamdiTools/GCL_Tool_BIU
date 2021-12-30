"""
GCL - Global Coordination Level Python Tool.
November 2021.
Code written by Omer Hamdi - omerhamdilf2@gmail.com.
"""
# Imports:
import gcl_library as gcl_lib

# check if input file/s exists:
import sys
from pathlib import Path
from os import path
# for threading:
import threading
# math tools:
import numpy as np
from scipy.spatial.distance import cdist
import math
import random
# gcl library import:
from matplotlib import pyplot as plt


def main_gcl_start(files_arr, jack_knifes=70, num_divisions=10, jack_knife_percentage=0.8, task='jackknife'):
    """
    The main gcl start function to generate and present the histograms.
    param files_arr: Array of .CSV files.
    param jack_knifes: Number of jackknives for calculation.
    param num_divisions: Number of random gene division for calculation.
    param jack_knife_percentage: Percentage of cells to choose for jackknife realization.
    param task: String, either 'jackknife' or 'regular_calc' for the requested task.
    return: none.
    """
    for file in files_arr:
        if not path.exists(file):
            raise NameError('Invalid File/s!')
    result_arr, file_names = [], []
    for file in files_arr:
        csv_mat = np.genfromtxt(file, delimiter=',')
        file_names.append(Path(file).stem)
        if task == 'jackknife':
            result_arr.append(gcl_lib.jackknife(csv_mat, jack_knifes, jack_knife_percentage, num_divisions))
        elif task == 'regular_calc':
            result_arr.append(gcl_lib.gcl(csv_mat, num_divisions))
            print('GCL value of ' + file_names[-1] + ' is: ' + str(result_arr[-1]))
        else:
            raise NameError('Invalid Task!')
    # plotting the histogram/s:
    plt.figure(figsize=(10, 10))
    if task == 'jackknife':
        for result in range(len(result_arr)):
            plt.hist(result_arr[result], 8, density=False, edgecolor='black', label=file_names[result], alpha=.8)
        plot_title = 'GCL ' + (
            'jackknife ' if task == 'jackknife' else 'Regular Calculation ') + 'Hist. with: ' + str(
            jack_knifes) + ' realizations' + ', ' + str(jack_knife_percentage * 100) + '%'
        plt.title(plot_title, fontsize=26)
    # plot settings:
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('GCL', fontsize=32)
    plt.ylabel('Iterations', fontsize=32)
    plt.rc('legend', fontsize=26)
    plt.legend(loc='upper right')
    plt.show()


if __name__ == '__main__':
    '''
    Run main_gcl_start with received arguments or default values.
    * Input: (as arguments, in the following order. for default write: 'default')
    *       1) data_arr [mandatory] - files array data with all the files in the current folder to make the GCL on.
    *       2) jackknife [optional] - Number of jackknives iterations for calculation.
    *       3) jackknife_percentage [optional] - Percentage of cells to choose for jackknife realization.
    *       4) task_option [optional] - either 'jackknife' or 'regular_calc' for the requested task.
    *       5) num_division [optional] - Number of random gene division iterations for calculation.
    '''

    if len(sys.argv) == 1:
        raise Exception("not enough arguments")
    data_arr = sys.argv[1].split(',')

    jack_knife = 70
    if len(sys.argv) > 2:
        jack_knife = 70 if sys.argv[2] == 'default' else int(sys.argv[2])

    jackknife_percentage = 0.8
    if len(sys.argv) > 3:
        jackknife_percentage = 0.8 if sys.argv[3] == 'default' else float(sys.argv[3])

    task_option = 'jackknife'
    if len(sys.argv) > 4:
        task_option = 'jackknife' if sys.argv[4] == 'default' else sys.argv[4]

    num_division = 10
    if len(sys.argv) > 5:
        num_division = 10 if sys.argv[5] == 'default' else int(sys.argv[5])

    main_gcl_start(data_arr, jack_knife, num_division, jackknife_percentage, task_option)
