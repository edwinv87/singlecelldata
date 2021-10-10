import pandas as pd
import numpy as np

import gzip
import os.path
import scipy.sparse
from scipy.sparse import csr_matrix

def ReadCSV(path, dset_name):

    data_path = path + dset_name + '/' + dset_name + "_data.csv"
    celldata_path = path + dset_name + '/' + dset_name + "_celldata.csv"
    genedata_path = path + dset_name + '/' + dset_name + "_genedata.csv"

    data = pd.read_csv(data_path, index_col=0)
    celldata = pd.read_csv(celldata_path, index_col=0)
    genedata = pd.read_csv(genedata_path, index_col = 0)

    return data, genedata, celldata





def Read10X(path):

    with open(path + '/matrix.mtx', 'r') as f:
        while True:
            header = f.readline()
            if not header.startswith('%'):
                break
        header = header.rstrip().split()
        n_genes, n_cells = int(header[0]), int(header[1])

        data, i, j = [], [], []
        for line in f:
            fields = line.rstrip().split()
            data.append(float(fields[2]))
            i.append(int(fields[1])-1)
            j.append(int(fields[0])-1)
        X = csr_matrix((data, (i, j)), shape=(n_cells, n_genes))

    genes = []
    with open(path + '/genes.tsv', 'r') as f:
        for line in f:
            fields = line.rstrip().split()
            genes.append(fields[1])
    assert(len(genes) == n_genes)

    barcodes = []
    with open(path + '/barcodes.tsv', 'r') as f:
        for line in f:
            fields = line.rstrip().split()
            barcodes.append(fields[1])
    assert(len(barcodes) == n_cells)

    return X, np.array(genes), np.array(barcodes)