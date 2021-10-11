# Imports:
import numpy as np
import math
import random
# alternative for 'pdist2' in matlab (euclidean distance):
from scipy.spatial.distance import cdist


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


def bcdcorr_calculation(data):
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


def gcl(data, num_divisions=100):
    """
    calculate the GCL of the data (num_divisions) times and return it as a list (vector).
    :param data: Input data such that [num_genes, num_cells] = size(data).
    :param num_divisions: Number of random gene division to calculate (default: 100)..
    :return: The GCL of the data - a number (nan mean value).
    """
    gcl_output = []
    for i in range(num_divisions):
        gcl_output.append(bcdcorr_calculation(data))
    return np.nanmean(gcl_output)


def bootstrap(data, boot_straps=100, num_divisions=100, choose_percentage=0.8):
    """
    bootstrapping the data - instead of repeating values, unique choose_percentage values from the data.
    :param data: Input data.
    :param boot_straps: Number of iterations to calculate (default: 100).
    :param num_divisions: Number of divisions for the gcl calculation (default: 100).
    :param choose_percentage: The percentage of the cells to calculate in each iteration (default: 0.8).
    :return:
    """
    boot_strap_arr = []
    for i in range(boot_straps):
        delete_arr = random.sample(range(0, len(data[0])), round(len(data[0]) * (1 - choose_percentage)))
        boot_strap_arr.append(gcl(np.delete(data, delete_arr, 1), num_divisions))
    return boot_strap_arr
