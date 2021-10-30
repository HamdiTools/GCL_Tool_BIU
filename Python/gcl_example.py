# Imports:
import gcl_library as gcl_lib

import sys
import numpy as np
# check if input file/s exists:
from pathlib import Path
from os import path
import threading
from scipy.spatial.distance import cdist
import math
import random
# gcl library import:
from matplotlib import pyplot as plt
import time


def main_gcl_start(files_arr, boot_straps=100, num_divisions=100, boot_strap_percentage=0.5, task='bootstrap'):
    '''
    The main gcl start function to generate and present the histograms.
    :param files_arr: Array of .CSV files.
    :param boot_straps: Number of bootstraps for calculation.
    :param num_divisions: Number of random gene division for calculation.
    :param boot_strap_percentage: Percentage of cells to choose for bootstrapping.
    :param task: String, either 'bootsrtap' or 'regular_calc' for the requested task.
    :return: none.
    '''
    for file in files_arr:
        if not path.exists(file):
            raise NameError('Invalid File/s!')
    result_arr, file_names = [], []
    for file in files_arr:
        csv_mat = np.genfromtxt(file, delimiter=',')
        file_names.append(Path(file).stem)
        if task == 'bootstrap':
            result_arr.append(gcl_lib.bootstrap(csv_mat, boot_straps, num_divisions, boot_strap_percentage))
        elif task == 'regular_calc':
            result_arr.append(gcl_lib.gcl(csv_mat, num_divisions))
            print('GCL value of ' + file_names[-1] + ' is: ' + str(result_arr[-1]))
        else:
            raise NameError('Invalid Task!')
    # plotting the histogram:
    if task == 'bootstrap':
        for result in range(len(result_arr)):
            plt.hist(result_arr[result], 8, density=False, edgecolor='black', label=file_names[result], alpha=.8)
        plot_title = 'GCL - ' + (
            'BootStrap ' if task == 'bootstrap' else 'Regular Calculation ') + 'Histogram with: ' + str(
            num_divisions) + ' iterations' + (', ' + str(boot_straps) + ' Boot Straps' + ', ' + str(
            boot_strap_percentage * 100) + '%' if task == 'bootstrap' else '')
        plt.title(plot_title)
        plt.xlabel('GCL')
        plt.ylabel('Iterations')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    '''
    run main_gcl_start with received arguments or default values.
    *Input: (as arguments, in the following order. for default write: 'default')
    *       1) data_arr [mandatory] - files array data with all the files in the current folder to make the GCL on.
    *       2) num_division [optional] - Number of random gene division iterations for calculation.
    *       3) bootsrtaps [optional] - Number of boot straps iterations for calculation.
    *       4) bootstrap_percentage [optional] - Percentage of cells to choose for bootstrapping.
    *       5) task_option [optional] - either 'bootsrtap' or 'regular_calc' for the requested task.
    '''
    for k in range(7):
        t = time.time()
        if len(sys.argv) == 1:
            raise Exception("not enough arguments")
        data_arr = sys.argv[1].split(',')
        num_division = 100
        if len(sys.argv) > 2:
            num_division = 100 if sys.argv[2] == 'default' else int(sys.argv[2])
        boot_strap = 100
        if len(sys.argv) > 3:
            boot_strap = 100 if sys.argv[3] == 'default' else int(sys.argv[3])
        bootstrap_percentage = 0.8
        if len(sys.argv) > 4:
            bootstrap_percentage = 0.8 if sys.argv[4] == 'default' else float(sys.argv[4])
        task_option = 'bootstrap'
        if len(sys.argv) > 5:
            task_option = 'bootstrap' if sys.argv[5] == 'default' else sys.argv[5]
        main_gcl_start(data_arr, boot_strap, num_division, bootstrap_percentage, task_option)
        print('Time for all proccess: ' + str(time.time() - t)+'seconds')
