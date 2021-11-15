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


def main_gcl_start(files_arr, boot_straps=70, num_divisions=10, boot_strap_percentage=0.8, task='bootstrap'):
    '''
    The main gcl start function to generate and present the histograms.
    :param files_arr: Array of .CSV files.
    :param boot_straps: Number of bootstraps for calculation.
    :param num_divisions: Number of random gene division for calculation.
    :param boot_strap_percentage: Percentage of cells to choose for bootstrapping.
    :param task: String, either 'bootstrap' or 'regular_calc' for the requested task.
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
            result_arr.append(gcl_lib.bootstrap(csv_mat, boot_straps, boot_strap_percentage, num_divisions))
        elif task == 'regular_calc':
            result_arr.append(gcl_lib.gcl(csv_mat, num_divisions))
            print('GCL value of ' + file_names[-1] + ' is: ' + str(result_arr[-1]))
        else:
            raise NameError('Invalid Task!')
    # plotting the histogram/s:
    if task == 'bootstrap':
        bin_width = 0.015
        overall_min, overall_max = 0.2, 0.550000000001
        for result in range(len(result_arr)):
            plt.hist(result_arr[result], density=False, edgecolor='black', label=file_names[result], alpha=.8,
                     bins=np.arange(min(result_arr[result]), max(result_arr[result]) + bin_width, bin_width))
        plot_title = 'GCL - ' + (
            'BootStrap ' if task == 'bootstrap' else 'Regular Calculation ') + 'Histogram with: ' + str(
            boot_straps) + ' Boot Straps' + ', ' + str(boot_strap_percentage * 100) + '%'
        # plt.title(plot_title)
        plt.xticks(np.arange(overall_min, overall_max, 0.025 * 4), fontsize=20)
        plt.yticks(np.arange(0, 50, 10),fontsize=20)
    plt.xlabel('GCL', fontsize=24)
    plt.ylabel('Iterations', fontsize=24)
    plt.rc('legend', fontsize='x-large')
    plt.legend(loc='upper right')
    plt.show()


if __name__ == '__main__':
    '''
    Run main_gcl_start with received arguments or default values.
    * Input: (as arguments, in the following order. for default write: 'default')
    *       1) data_arr [mandatory] - files array data with all the files in the current folder to make the GCL on.
    *       2) bootstraps [optional] - Number of boot straps iterations for calculation.
    *       3) bootstrap_percentage [optional] - Percentage of cells to choose for bootstrapping.
    *       4) task_option [optional] - either 'bootstrap' or 'regular_calc' for the requested task.
    *       5) num_division [optional] - Number of random gene division iterations for calculation.
    '''

    if len(sys.argv) == 1:
        raise Exception("not enough arguments")
    data_arr = sys.argv[1].split(',')

    boot_strap = 70
    if len(sys.argv) > 2:
        boot_strap = 70 if sys.argv[2] == 'default' else int(sys.argv[2])

    bootstrap_percentage = 0.8
    if len(sys.argv) > 3:
        bootstrap_percentage = 0.8 if sys.argv[3] == 'default' else float(sys.argv[3])

    task_option = 'bootstrap'
    if len(sys.argv) > 4:
        task_option = 'bootstrap' if sys.argv[4] == 'default' else sys.argv[4]

    num_division = 10
    if len(sys.argv) > 5:
        num_division = 10 if sys.argv[5] == 'default' else int(sys.argv[5])

    main_gcl_start(data_arr, boot_strap, num_division, bootstrap_percentage, task_option)
