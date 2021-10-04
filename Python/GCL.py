'''
Calculate the GCL of the input data with bootstrap
Guy Amit, guy1.amit@gmail.com, Orr Levy, Dana Vaknin, Tom Snir, Sol
Efroni, Peter Castaldi, Yang-Yu Liu, Haim Cohen, Amir Bashan.
Based on Bias Corrected Distance Correlation Szekely, G. J., & Rizzo,
M. L. (2013). The distance correlation t-test of independence in
high dimension. Journal of Multivariate Analysis, 117, 193-213.
data - Input data such that [num_genes, num_cells] = size(data)
num_division - Number of random gene division for calculation.
gcl_output - The GCL of the data
'''

# Imports:
import numpy as np
import matplotlib.pyplot as plt
import math
import random
# for receiving arguments:
import sys
# alternative for 'pdist2' in matlab (euclidean distance):
from scipy.spatial.distance import cdist
# check if input file/s exists:
from pathlib import Path
from os import path


def vn(Aij, Bij, cells):
    """
    Calculate the variance.
    :param Aij: first matrix.
    :param Bij: second matrix.
    :param cells: number of cells.
    :return: new Aij matrix.
    """
    return np.dot((1 / (cells * (cells - 3))),
                  sum(sum(Aij * Bij)) - np.matmul(np.dot(cells / (cells - 2), np.diag(Aij).T), np.diag(Bij)))


def rn(Aij, Bij, cells):
    """
    calculate the covariance.
    :param Aij: first matrix.
    :param Bij: second matrix.
    :param cells: number of cells.
    :return: new Aij matrix.
    """
    return np.divide(vn(Aij, Bij, cells), math.sqrt((vn(Aij, Aij, cells) * vn(Bij, Bij, cells))))


def get_matrix(genes_from_data, cells):
    """
    calculate the euclidean distance of the genes_from_data ->d.
    calculate the mean and vector_mean of d.
    calculate using matrix operations the Aij matrix for the GCL calculation.
    :param genes_from_data: the specific genes from the data to calculate the matrix.
    :param cells: number of cells of each gene.
    :return: Aij matrix.
    """
    d = cdist(genes_from_data, genes_from_data, metric='euclidean')
    m, vector_m = np.mean(d), np.mean(d, axis=0).reshape(1, cells)
    # operations on the matrices:
    Aij = d - np.matmul(vector_m.T, np.ones((1, cells))) - np.dot(np.ones((cells, 1)), vector_m) + m - vector_m / cells
    np.fill_diagonal(Aij, vector_m - m)
    Aij = (cells / (cells - 1)) * Aij
    return Aij


def gcl_single_calculation(data):
    """
    split the data randomly to two equal parts, get their matrices and return their GCL.
    :param data: Input data.
    :return: GCL of the data.
    """
    num_genes, cells = len(data), len(data[0])
    random_genes = np.random.permutation(num_genes)
    first_half, second_half = random_genes[:math.floor(num_genes / 2)], random_genes[math.floor(num_genes / 2):]
    Aij1, Aij2 = get_matrix(np.transpose(data[first_half]), cells), get_matrix(np.transpose(data[second_half]), cells)
    # Calculating bcdcorr:
    return rn(Aij2, Aij1, cells)


def gcl_calculator(data, num_divisions):
    """
    calculate the GCL of the data (num_divisions) times and return it as a list (vector).
    :param data: Input data such that [num_genes, num_cells] = size(data).
    :param num_divisions: Number of random gene division to calculate.
    :return: The GCL of the data.
    """
    gcl_output = []
    for i in range(num_divisions):
        gcl_output.append(gcl_single_calculation(data))
    return gcl_output


def bootstrap(data, boot_straps, choose_percentage):
    """
    bootstrapping the data - instead of repeating values, unique choose_percentage values from the data.
    :param data: Input data.
    :param boot_straps: Number of iterations to calculate.
    :param choose_percentage: The percentage of the cells to calculate in each iteration.
    :return:
    """
    boot_strap_arr = []
    for i in range(boot_straps):
        delete_arr = random.sample(range(0, len(data[0])), round(len(data[0]) * (1 - choose_percentage)))
        boot_strap_arr.append(gcl_single_calculation(np.delete(data, delete_arr, 1)))
    return boot_strap_arr


def main_gcl_start(files_arr, num_divisions=100, boot_strap_percentage=0.5, task='bootstrap'):
    '''
    The main gcl start function to generate and present the histograms.
    :param files_arr: Array of .CSV files.
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
            result_arr.append(bootstrap(csv_mat, num_divisions, boot_strap_percentage))
        elif task == 'regular_calc':
            result_arr.append(gcl_calculator(csv_mat, num_divisions))
        else:
            raise NameError('Invalid Task!')
    for result in range(len(result_arr)):
        plt.hist(result_arr[result], 8, density=False, edgecolor='black', label=file_names[result], alpha=.8)
    plot_title = 'GCL - ' + ('BootStrap ' if task == 'bootstrap' else 'Regular Calculation ') + 'Histogram with ' + str(
        num_divisions) + ' iterations' + (', ' + str(boot_strap_percentage) + '%' if task == 'bootstrap' else '')
    plt.title(plot_title)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    '''
    run main_gcl_start with received arguments or default values.
    *Input: (as arguments, in the following order. for default write: 'default')
    *       1) data_arr [mandatory] - files array data with all the files in the current folder to make the GCL on.
    *       2) num_division [optional] - Number of random gene division/ bootsrtap iterations for calculation.
    *       3) bootstrap_percentage [optional] - Percentage of cells to choose for bootstrapping.
    *       4) task_option [optional] - either 'bootsrtap' or 'regular_calc' for the requested task.
    '''
    data_arr = sys.argv[1].split(',')
    num_division = 100
    if len(sys.argv) > 2:
        num_division = 100 if sys.argv[2] == 'default' else int(sys.argv[2])
    bootstrap_percentage = 0.8
    if len(sys.argv) > 3:
        bootstrap_percentage = 0.8 if sys.argv[3] == 'default' else float(sys.argv[3])
    task_option = 'bootstrap'
    if len(sys.argv) > 4:
        task_option = 'bootstrap' if sys.argv[4] == 'default' else sys.argv[4]
    main_gcl_start(data_arr, num_division, bootstrap_percentage, task_option)
